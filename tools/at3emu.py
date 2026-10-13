"""Run psp_at3tool.exe (32-bit Windows PE) on any host via the Unicorn CPU emulator.

Only ~25 MSVCR80 imports are used, so they are shimmed in Python (file I/O, heap,
printf) or as tiny native x87 stubs (floor/pow/asin/ftol) to keep FPU-stack semantics
exact. The FPU control word is set to 0x27F, matching the MSVC CRT start-up call
_controlfp_s(NULL, _PC_53, _MCW_PC).

CLI:   python at3emu.py -e in.wav out.at3 -br 132
API:   emu = Emu(); emu.run_main(['-e', ...]);  emu.call(addr, *args); emu.hook(addr, fn)
"""
import os, sys, struct, math
import pefile
from unicorn import Uc, UC_ARCH_X86, UC_MODE_32, UC_HOOK_CODE, UC_PROT_ALL
from unicorn.x86_const import *

EXE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'private', 'psp_at3tool.exe')

STUB_BASE = 0x70000000      # import stubs
SCRATCH = 0x70010000        # scratch for float returns
IOB_BASE = 0x70020000       # fake FILE structs (__iob_func), 0x20 bytes each
FILE_BASE = 0x70021000      # fake FILE structs for fopen_s
STRS = 0x70030000           # argv strings
STACK_TOP = 0x7f000000
STACK_SIZE = 16 << 20
HEAP_BASE = 0x10000000
HEAP_SIZE = 512 << 20
RET_MAGIC = 0x70000ff0      # return address that ends an emu.call()


class ExitProgram(Exception):
    pass


def _asm_float_stub(kind):
    """Native x87 stubs; the 'nop' (0x90) is where the Python hook computes the value."""
    if kind == 'pow':      # ST0 = y, ST1 = x  -> pow(x, y)
        return bytes([0x83, 0xEC, 0x10, 0xDD, 0x5C, 0x24, 0x08, 0xDD, 0x1C, 0x24, 0x90,
                       0xDD, 0x04, 0x24, 0x83, 0xC4, 0x10, 0xC3]), 10
    if kind == 'asin':     # ST0 = x -> asin(x)
        return bytes([0x83, 0xEC, 0x08, 0xDD, 0x1C, 0x24, 0x90,
                       0xDD, 0x04, 0x24, 0x83, 0xC4, 0x08, 0xC3]), 6
    if kind == 'floor':    # cdecl double arg -> ST0
        return bytes([0x90, 0xDD, 0x05]) + struct.pack('<I', SCRATCH) + bytes([0xC3]), 0
    if kind == 'ftol':     # truncate ST0 -> EDX:EAX, pop
        return bytes([0x83, 0xEC, 0x0C, 0xD9, 0x7C, 0x24, 0x08, 0x66, 0x8B, 0x44, 0x24, 0x08,
                       0x80, 0xCC, 0x0C, 0x66, 0x89, 0x44, 0x24, 0x0A, 0xD9, 0x6C, 0x24, 0x0A,
                       0xDF, 0x3C, 0x24, 0xD9, 0x6C, 0x24, 0x08, 0x8B, 0x04, 0x24, 0x8B, 0x54,
                       0x24, 0x04, 0x83, 0xC4, 0x0C, 0xC3]), None
    raise ValueError(kind)


class Emu:
    def __init__(self, exe=EXE, quiet=False):
        self.quiet = quiet
        self.out = []                      # captured stdout text
        self.files = {}                    # FILE* -> python file object
        self.eof = set()                   # id(file) with the C EOF indicator set
        self.next_file = FILE_BASE
        self.heap_ptr = HEAP_BASE
        self.allocs = {}
        self.hooks = {}
        mu = self.mu = Uc(UC_ARCH_X86, UC_MODE_32)
        pe = self.pe = pefile.PE(exe)
        base = pe.OPTIONAL_HEADER.ImageBase
        mu.mem_map(base, (pe.OPTIONAL_HEADER.SizeOfImage + 0xfff) & ~0xfff, UC_PROT_ALL)
        mu.mem_write(base, pe.header[:0x1000])
        for s in pe.sections:
            data = s.get_data()
            mu.mem_write(base + s.VirtualAddress, data[:max(s.Misc_VirtualSize, len(data))])
        mu.mem_map(STUB_BASE, 0x40000, UC_PROT_ALL)
        mu.mem_map(STACK_TOP - STACK_SIZE, STACK_SIZE, UC_PROT_ALL)
        mu.mem_map(HEAP_BASE, HEAP_SIZE, UC_PROT_ALL)
        self._install_imports()
        mu.hook_add(UC_HOOK_CODE, self._on_code, begin=STUB_BASE, end=STUB_BASE + 0xffff)

    # ---------------------------------------------------------------- memory helpers
    def r32(self, a): return struct.unpack('<I', self.mu.mem_read(a, 4))[0]
    def w32(self, a, v): self.mu.mem_write(a, struct.pack('<I', v & 0xffffffff))
    def rf64(self, a): return struct.unpack('<d', self.mu.mem_read(a, 8))[0]
    def wf64(self, a, v): self.mu.mem_write(a, struct.pack('<d', v))
    def read(self, a, n): return bytes(self.mu.mem_read(a, n))
    def write(self, a, b): self.mu.mem_write(a, bytes(b))

    def cstr(self, a):
        out = bytearray()
        while True:
            c = self.mu.mem_read(a, 1)[0]
            if c == 0: return out.decode('latin-1')
            out.append(c); a += 1

    def arg(self, i):  # cdecl arg i at hook time (return address at [esp])
        return self.r32(self.mu.reg_read(UC_X86_REG_ESP) + 4 + 4 * i)

    def malloc(self, n):
        n = max(n, 1)
        p = (self.heap_ptr + 15) & ~15
        self.heap_ptr = p + n
        if self.heap_ptr > HEAP_BASE + HEAP_SIZE: raise MemoryError('emulated heap exhausted')
        self.allocs[p] = n
        return p

    # ---------------------------------------------------------------- imports
    def _install_imports(self):
        self.stubs = {}
        cur = STUB_BASE
        for entry in self.pe.DIRECTORY_ENTRY_IMPORT:
            for imp in entry.imports:
                name = imp.name.decode() if imp.name else 'ord%d' % imp.ordinal
                kind = {'_CIpow': 'pow', '_CIasin': 'asin', 'floor': 'floor', '_ftol': 'ftol'}.get(name)
                if kind:
                    code, hook_off = _asm_float_stub(kind)
                    self.mu.mem_write(cur, code)
                    if hook_off is not None: self.stubs[cur + hook_off] = name
                else:
                    self.mu.mem_write(cur, b'\x90\xc3')  # nop (hooked) ; ret
                    self.stubs[cur] = name
                self.w32(imp.address, cur)
                cur += 0x40
        self.mu.mem_write(RET_MAGIC, b'\x90')
        self.stubs[RET_MAGIC] = '__return__'
        self.w32(IOB_BASE + 0x0, 0)

    def _on_code(self, mu, addr, size, _):
        name = self.stubs.get(addr)
        if name is None: return
        if name == '__return__':
            mu.emu_stop(); return
        fn = getattr(self, 'imp_' + name, None)
        if fn is None:
            raise NotImplementedError('import %s not shimmed' % name)
        ret = fn()
        if ret is not None:
            mu.reg_write(UC_X86_REG_EAX, ret & 0xffffffff)

    # float stubs: operands spilled to [esp] by the native prologue
    def imp__CIpow(self):
        esp = self.mu.reg_read(UC_X86_REG_ESP)
        x, y = self.rf64(esp), self.rf64(esp + 8)
        try: r = math.pow(x, y)
        except (OverflowError, ValueError): r = float('nan') if x < 0 else float('inf')
        self.wf64(esp, r)

    def imp__CIasin(self):
        esp = self.mu.reg_read(UC_X86_REG_ESP)
        self.wf64(esp, math.asin(self.rf64(esp)))

    def imp_floor(self):
        esp = self.mu.reg_read(UC_X86_REG_ESP)
        self.wf64(SCRATCH, math.floor(self.rf64(esp + 4)))

    # heap / memory
    def imp_malloc(self): return self.malloc(self.arg(0))
    def imp_calloc(self):
        n = self.arg(0) * self.arg(1); p = self.malloc(n)
        self.mu.mem_write(p, b'\0' * n); return p
    def imp_free(self): return 0
    def imp_memset(self):
        d, c, n = self.arg(0), self.arg(1) & 0xff, self.arg(2)
        if n: self.mu.mem_write(d, bytes([c]) * n)
        return d
    def imp_memcpy(self):
        d, s, n = self.arg(0), self.arg(1), self.arg(2)
        if n: self.mu.mem_write(d, bytes(self.mu.mem_read(s, n)))
        return d
    imp_memmove = imp_memcpy

    # stdio
    def _fmt(self, fmt, argi):
        out, i = [], 0
        while i < len(fmt):
            c = fmt[i]
            if c != '%':
                out.append(c); i += 1; continue
            j = i + 1
            while fmt[j] in '-+ #0123456789.lh': j += 1
            spec, conv = fmt[i:j], fmt[j]
            spec = spec.replace('l', '').replace('h', '')
            if conv == '%': out.append('%')
            elif conv == 's': out.append((spec + 's') % self.cstr(self.arg(argi))); argi += 1
            elif conv in 'di':
                v = self.arg(argi); v = v - (1 << 32) if v & 0x80000000 else v
                out.append((spec + 'd') % v); argi += 1
            elif conv in 'uxXc':
                out.append((spec + conv) % self.arg(argi)); argi += 1
            else:
                raise NotImplementedError('printf %' + conv)
            i = j + 1
        return ''.join(out)

    def _emit(self, s, err=False):
        self.out.append(s)
        if not self.quiet: (sys.stderr if err else sys.stdout).write(s)

    def imp_printf(self): self._emit(self._fmt(self.cstr(self.arg(0)), 1)); return 0
    def imp_fprintf(self): self._emit(self._fmt(self.cstr(self.arg(1)), 2), err=True); return 0
    def imp___iob_func(self): return IOB_BASE

    def imp_fopen_s(self):
        pf, name, mode = self.arg(0), self.cstr(self.arg(1)), self.cstr(self.arg(2))
        try: f = open(name, mode)
        except OSError: return 2
        h = self.next_file; self.next_file += 0x20
        self.files[h] = f; self.w32(pf, h); return 0
    def imp_fclose(self):
        f = self.files.pop(self.arg(0), None)
        if f: f.close()
        return 0
    def imp_fread(self):
        buf, sz, cnt, f = self.arg(0), self.arg(1), self.arg(2), self.files[self.arg(3)]
        if sz == 0 or cnt == 0: return 0
        data = f.read(sz * cnt)
        if len(data) < sz * cnt: self.eof.add(id(f))
        if data: self.mu.mem_write(buf, data)
        return len(data) // sz
    def imp_fwrite(self):
        buf, sz, cnt, f = self.arg(0), self.arg(1), self.arg(2), self.files[self.arg(3)]
        if sz * cnt: f.write(bytes(self.mu.mem_read(buf, sz * cnt)))
        return cnt
    def imp_fseek(self):
        f, off, whence = self.files[self.arg(0)], self.arg(1), self.arg(2)
        off = off - (1 << 32) if off & 0x80000000 else off
        f.seek(off, whence); self.eof.discard(id(f)); return 0
    def imp_ftell(self): return self.files[self.arg(0)].tell()
    def imp_feof(self): return 1 if id(self.files[self.arg(0)]) in self.eof else 0
    def imp_fgetc(self):
        f = self.files[self.arg(0)]
        b = f.read(1)
        if not b:
            self.eof.add(id(f))
            return 0xffffffff
        return b[0]
    def imp_atoi(self):
        s = self.cstr(self.arg(0)).strip(); n = ''
        for k, ch in enumerate(s):
            if ch.isdigit() or (k == 0 and ch in '+-'): n += ch
            else: break
        try: return int(n)
        except ValueError: return 0
    def imp_exit(self):
        self.exit_code = self.arg(0); self.mu.emu_stop(); self._exited = True

    # ---------------------------------------------------------------- execution
    def call(self, addr, *args):
        """cdecl call into the image; returns EAX."""
        mu = self.mu
        sp = getattr(self, '_sp', STACK_TOP - 0x1000)
        sp -= 4 * len(args)
        for k, a in enumerate(args): self.w32(sp + 4 * k, a)
        sp -= 4; self.w32(sp, RET_MAGIC)
        mu.reg_write(UC_X86_REG_ESP, sp)
        mu.reg_write(UC_X86_REG_FPCW, 0x27F)
        self._exited = False
        mu.emu_start(addr, RET_MAGIC + 1)
        return mu.reg_read(UC_X86_REG_EAX)

    def hook(self, addr, fn):
        """Call fn(emu) whenever execution reaches addr."""
        self.mu.hook_add(UC_HOOK_CODE, lambda mu, a, s, u: fn(self), begin=addr, end=addr)

    def put_str(self, s):
        p = self.malloc(len(s) + 1); self.mu.mem_write(p, s.encode('latin-1') + b'\0'); return p

    def run_main(self, argv):
        argv = ['at3tool'] + list(argv)
        ptrs = [self.put_str(a) for a in argv]
        arr = self.malloc(4 * (len(ptrs) + 1))
        for k, p in enumerate(ptrs + [0]): self.w32(arr + 4 * k, p)
        rc = self.call(0x4012f0, len(argv), arr)
        for f in self.files.values(): f.close()
        self.files.clear()
        return getattr(self, 'exit_code', rc) if self._exited else rc


if __name__ == '__main__':
    e = Emu()
    sys.exit(e.run_main(sys.argv[1:]))
