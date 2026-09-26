# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""Process priority for the hsp runs (2026-09-25; Moki&Julio).

geomd.lower_priority() (and the cube probe's pattern it copies) calls kernel32.SetPriorityClass(GetCurrentProcess(), ...)
through ctypes WITHOUT restype / argtypes. On 64-bit Python the pseudo-handle -1 returned by GetCurrentProcess is cut to
a 32-bit int, so SetPriorityClass receives a bad handle and returns 0: the call is a NO-OP (measured 2026-09-25:
SetPriorityClass -> 0, GetPriorityClass -> 0). lower() does it properly: psutil if installed, else ctypes with the
Win32 signatures declared. It returns the priority class read back from the OS (BELOW_NORMAL = 0x4000 = 16384).
"""
import os


def lower():
    try:
        import psutil
        p = psutil.Process(os.getpid()); p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        return int(p.nice())
    except Exception:
        pass
    try:
        import ctypes
        from ctypes import wintypes as wt
        k = ctypes.WinDLL('kernel32', use_last_error=True)
        k.GetCurrentProcess.restype = wt.HANDLE; k.GetCurrentProcess.argtypes = []
        k.SetPriorityClass.argtypes = [wt.HANDLE, wt.DWORD]; k.SetPriorityClass.restype = wt.BOOL
        k.GetPriorityClass.argtypes = [wt.HANDLE]; k.GetPriorityClass.restype = wt.DWORD
        h = k.GetCurrentProcess(); k.SetPriorityClass(h, 0x4000)
        return int(k.GetPriorityClass(h))
    except Exception:
        return -1


def patch_geomd():
    """Make geomd.lower_priority (called by the reused search / seeds code) the working version."""
    import geomd
    geomd.lower_priority = lower


if __name__ == '__main__':
    print('priority class after lower():', lower())
