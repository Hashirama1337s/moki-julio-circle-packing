# Moki&Julio circle-packing records - https://github.com/Hashirama1337s/moki-julio-circle-packing
"""BelowNormal priority for every chx process (same method as the ball solver's prio.py: psutil, else ctypes with declared Win32
signatures).  lower() returns the priority class READ BACK from the OS (BELOW_NORMAL = 0x4000 = 16384).  Moki&Julio."""
import os

BELOW_NORMAL = 0x4000


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
        h = k.GetCurrentProcess(); k.SetPriorityClass(h, BELOW_NORMAL)
        return int(k.GetPriorityClass(h))
    except Exception:
        return -1


if __name__ == '__main__':
    print('priority class after lower():', lower())
