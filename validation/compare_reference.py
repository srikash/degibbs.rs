#!/usr/bin/env python3
"""Compare degibbs with a compiled RPG Python module (the DESIGNER-v2 `rpg*.so`) on the Shepp-Logan phantom.

  python3 validation/compare_reference.py --exe target/release/degibbs --rpg-path /dir/with/rpg.so --out /tmp/cmp

Prints the relative RMS difference for every acquired and correction partial-Fourier factor. Expect 2e-5 or
better for 1, 7/8 and 6/8. For 5/8 magnitude input the DESIGNER-v2 module skips the sub-image step, so it differs.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phantom import PFS, Setup, run


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--exe', required=True); p.add_argument('--rpg-path', required=True); p.add_argument('--out', required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sys.path.insert(0, a.rpg_path)
    import rpg
    n = 96
    st = Setup(n, 8, 0, 'smooth')
    for acq in PFS:
        img = st.acquire(acq, 'hi', None, None)
        for cp in PFS:
            arr = np.ascontiguousarray(np.broadcast_to(img, (1, 1, n, n))).astype(np.float64).copy()
            ref = np.asarray(rpg.unring(arr, minW=1, maxW=3, nsh=20, pfv=cp, pfdimf=True, phase_flag=False)[0]).reshape(n, n)
            mine = run(a.exe, a.out, img, cp, 'cmp')
            rel = np.sqrt(np.mean((mine - ref) ** 2)) / np.sqrt(np.mean(img ** 2))
            print(f'acquired {acq:.3f}  correction {cp:.3f}  relative RMS difference {rel:.2e}')


if __name__ == '__main__':
    main()
