#!/usr/bin/env python3
"""Partial-Fourier validation on a Shepp-Logan phantom.

Truncates the k-space of a high-resolution phantom (optionally with a smooth phase and noise), zeroes one end of
the phase-encode (y) lines to mimic partial Fourier, takes the magnitude, runs `degibbs` with every --pf setting,
and reports errors against the voxel-averaged truth. Averages over sub-voxel offsets of the object.

  python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --figures docs

Needs numpy, scipy, nibabel and matplotlib (for --figures).
"""
import argparse, json, os, subprocess, sys
import numpy as np, nibabel as nib
from scipy import ndimage as ndi

# Modified Shepp-Logan (Toft 1996): amplitude, semi-axes a, b, centre x0, y0, angle in degrees; coordinates in [-1, 1].
SL = [(1, .69, .92, 0, 0, 0), (-.8, .6624, .874, 0, -.0184, 0), (-.2, .11, .31, .22, 0, -18), (-.2, .16, .41, -.22, 0, 18),
      (.1, .21, .25, 0, .35, 0), (.1, .046, .046, 0, .1, 0), (.1, .046, .046, 0, -.1, 0), (.1, .046, .023, -.08, -.605, 0),
      (.1, .023, .023, 0, -.606, 0), (.1, .023, .046, .06, -.605, 0)]
PFS = (1, 7 / 8, 6 / 8, 5 / 8)
NAMES = ['1', '7/8', '6/8', '5/8']
BRAIN = 20.0  # brain intensity once the skull is scaled to 100


def coords(h):
    y, x = (np.mgrid[0:h, 0:h] + .5) / h * 2 - 1
    return -y, x  # rows run downwards, y points up


def ellipse(h, k):
    y, x = coords(h)
    _, a, b, x0, y0, t = SL[k]
    t = np.deg2rad(t); c, s = np.cos(t), np.sin(t)
    u = (x - x0) * c + (y - y0) * s; v = -(x - x0) * s + (y - y0) * c
    return (u / a) ** 2 + (v / b) ** 2 <= 1


def shepp_logan(h):
    return 100 * sum(SL[k][0] * ellipse(h, k) for k in range(len(SL)))


class Setup:
    def __init__(self, n, r, offset, phase):
        self.n, self.r, h = n, r, n * r
        self.h = h
        obj = shepp_logan(h)
        inner, outer = ellipse(h, 1), ellipse(h, 0)
        yy, xx = np.mgrid[0:h, 0:h] / h
        if phase == 'smooth':
            ph = np.exp(1j * (1.2 * np.sin(2 * np.pi * xx * .7) + np.cos(2 * np.pi * yy * .5)))
        elif phase == 'fast':   # phase turns over about 10 voxels, short against the 4-voxel 5/8 ringing interval
            ph = np.exp(1j * (1.5 * np.sin(2 * np.pi * xx * n / 10) + 1.5 * np.cos(2 * np.pi * yy * n / 10)))
        else:
            ph = 1.0
        sh = lambda a: np.roll(a, offset, axis=(0, 1))
        obj, inner, outer = sh(obj), sh(inner), sh(outer)
        ph = sh(ph) if phase != 'none' else ph
        k = np.fft.fftshift(np.fft.fft2(obj * ph)); c = h // 2
        self.k = k[c - n // 2:c + n // 2, c - n // 2:c + n // 2] / r ** 2
        # truth blocks are centred on the sample positions of the truncated image
        blk = lambda a: np.roll(a.astype(float), r // 2, axis=(0, 1)).reshape(n, r, n, r).mean((1, 3))
        self.truth = blk(obj)
        gy, gx = np.gradient(self.truth)
        edge = np.hypot(gx, gy) > 1e-6
        dist = ndi.distance_transform_edt(~ndi.binary_dilation(edge, iterations=1))
        far, near = dist > 4, dist <= 2
        in_in = blk(inner) > .999
        in_core = ndi.binary_erosion(in_in, iterations=3)       # inside the inner skull, away from it
        out_out = blk(outer) < .001
        skull = ndi.binary_dilation((blk(outer) > .001) & ~in_in, iterations=2) & ~in_core
        self.mask = {
            'far_brain': far & in_core & (self.truth > BRAIN / 2),   # brain tissue, no ventricles
            'far_background': far & out_out,
            'far_all': far,
            'near_brain': near & in_core,
            'skull': skull,
            'whole_brain': in_core,
            'whole_all': np.ones_like(far),
        }

    def acquire(self, pf, side, snr, rng):
        n = self.n
        m = np.ones(n, bool)
        keep = int(round(pf * n))
        if pf < 1:
            if side == 'hi':
                m[keep:] = False
            else:
                m[:] = True; m[:n - keep] = False
        k = self.k * m[:, None]
        if snr:
            sk = BRAIN / snr * n   # image-domain sigma per channel (real and imaginary) = brain / snr
            k = k + m[:, None] * (rng.normal(size=k.shape) + 1j * rng.normal(size=k.shape)) * sk
        return np.abs(np.fft.ifft2(np.fft.ifftshift(k)))


def run(exe, out, img, pf, tag):
    inp, res = f'{out}/{tag}.nii.gz', f'{out}/{tag}_c.nii.gz'
    nib.save(nib.Nifti1Image(np.repeat(img.T[:, :, None], 4, 2).astype(np.float32), np.eye(4)), inp)
    subprocess.run([exe, inp, res, '--pe-axis', 'j', '--pf', str(pf), '--no-json'], check=True, capture_output=True)
    return np.asanyarray(nib.load(res).dataobj)[:, :, 0].T.astype(float)


def evaluate(a, exe, out, offsets, side, phase, snr, seeds, n, r):
    metrics = ['far_brain', 'far_background', 'far_all', 'near_brain', 'skull', 'whole_brain', 'whole_all']
    # res[metric][acq, 0 = none / 1..4 = correction, offset*seed]
    res = {m: [] for m in metrics}
    for oi, off in enumerate(offsets):
        st = Setup(n, r, off, phase)
        for seed in seeds:
            rng = np.random.default_rng(1000 * oi + seed)
            cur = {m: np.zeros((4, 5)) for m in metrics}
            for i, pf in enumerate(PFS):
                img = st.acquire(pf, side, snr, rng)
                outs = [img] + [run(exe, out, img, cp, f'w_{i}') for cp in PFS]
                for j, o in enumerate(outs):
                    for m in metrics:
                        cur[m][i, j] = np.sqrt(np.mean((o - st.truth)[st.mask[m]] ** 2))
            for m in metrics:
                res[m].append(cur[m])
    return {m: np.array(v) for m, v in res.items()}


def show(res, title, counts):
    print(f'\n== {title} ==')
    print('pixels per mask (offset 0):', {m: int(c) for m, c in counts.items()})
    for m, v in res.items():
        mean = v.mean(0)
        print(f'-- {m}: mean over {len(v)} runs (rows: acquired; cols: none, --pf 1, 7/8, 6/8, 5/8)')
        for i in range(4):
            print(f'   {NAMES[i]:>4}', ' '.join(f'{x:7.3f}' for x in mean[i]))


def figures(args, st, exe, out, summary):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    n, truth, mask = args.n, st.truth, st.mask
    rng = np.random.default_rng(0)
    cases = []
    for i, pf in enumerate(PFS):
        img = st.acquire(pf, 'hi', None, rng)
        cases.append((img, run(exe, out, img, pf, f'qc{i}'), run(exe, out, img, 1, f'qc{i}f')))
    rm = lambda a, m: np.sqrt(np.mean((a - truth)[mask[m]] ** 2))
    fig = plt.figure(figsize=(17, 9.4), constrained_layout=True)
    gs = fig.add_gridspec(4, 9, width_ratios=[1, 1, 1, 1, .08, 1, 1, 1, .08])
    for i, (img, mt, fu) in enumerate(cases):
        for j, (a, t) in enumerate([(truth, 'truth'), (img, 'uncorrected'), (mt, f'corr {NAMES[i]}'), (fu, 'corr full-Fourier')]):
            ax = fig.add_subplot(gs[i, j]); ax.imshow(a, cmap='gray', vmin=0, vmax=truth.max()); ax.axis('off')
            ax.set_title(f'acq PF {NAMES[i]}: {t}', fontsize=9)
        for j, (a, t) in enumerate([(img, 'none'), (mt, 'matching'), (fu, 'full-F')]):
            ax = fig.add_subplot(gs[i, 5 + j]); ei = ax.imshow(np.abs(a - truth), cmap='magma', vmin=0, vmax=6); ax.axis('off')
            ax.set_title(f'|err| {t}\nfar brain {rm(a, "far_brain"):.2f}  far background {rm(a, "far_background"):.2f}', fontsize=9)
    fig.colorbar(ei, cax=fig.add_subplot(gs[:, 8]), label='|image - truth|  (clipped at 6)')
    fig.savefig(f'{args.figures}/qc_images.png', dpi=80)
    fig, axs = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    xc = n // 2 + 10
    for i, (img, mt, fu) in enumerate(cases):
        ax = axs.flat[i]
        ax.plot(truth[:, xc], 'k', lw=2, label='truth'); ax.plot(img[:, xc], 'r', lw=1, label='uncorrected')
        ax.plot(mt[:, xc], 'b', lw=1, label=f'corr {NAMES[i]}'); ax.plot(fu[:, xc], 'g', lw=1, label='corr full')
        ax.set_ylim(12, 28); ax.set_xlim(14, 82); ax.legend(fontsize=8); ax.set_xlabel('y (PE)')
        ax.set_title(f'acq PF {NAMES[i]}: PE profile, column x={xc}')
    for slot, (m, ttl) in zip((4, 5), (('far_brain', 'far-field RMSE, brain'), ('far_background', 'far-field RMSE, background'))):
        ax = axs.flat[slot]; M = summary[m][:, 1:]
        ax.imshow(M, cmap='viridis'); ax.set_xticks(range(4), NAMES); ax.set_yticks(range(4), NAMES)
        ax.set_xlabel('correction PF'); ax.set_ylabel('acquired PF'); ax.set_title(ttl)
        for i in range(4):
            for j in range(4):
                ax.text(j, i, f'{M[i, j]:.2f}', ha='center', va='center', color='w', fontsize=10)
    fig.savefig(f'{args.figures}/qc_profiles.png', dpi=80)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--exe', required=True, help='path to the degibbs binary')
    p.add_argument('--out', required=True, help='scratch directory for NIfTI files')
    p.add_argument('--n', type=int, default=96); p.add_argument('--r', type=int, default=8)
    p.add_argument('--offsets', type=int, default=8, help='number of sub-voxel object offsets to average')
    p.add_argument('--phase', choices=['smooth', 'fast', 'none'], default='smooth', help='phase map: slow, fast (about 10 voxels) or none')
    p.add_argument('--side', choices=['hi', 'lo'], default='hi', help='which end of k-space is left unsampled')
    p.add_argument('--snr', type=float, default=None, help='brain intensity divided by the noise sigma of each complex channel')
    p.add_argument('--seeds', type=int, default=1)
    p.add_argument('--figures', help='directory to write qc_images.png and qc_profiles.png (offset 0, noise-free)')
    p.add_argument('--json', help='write the mean tables to this file')
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.figures:
        os.makedirs(a.figures, exist_ok=True)
    offsets = [(k, (3 * k) % a.r) for k in range(a.offsets)] if a.offsets > 1 else [(0, 0)]
    res = evaluate(a, a.exe, a.out, offsets, a.side, a.phase, a.snr, range(a.seeds), a.n, a.r)
    counts = {m: v.sum() for m, v in Setup(a.n, a.r, 0, a.phase).mask.items()}
    show(res, f'phase={a.phase} side={a.side} snr={a.snr} offsets={len(offsets)} seeds={a.seeds}', counts)
    summary = {m: v.mean(0) for m, v in res.items()}
    # spread of the close comparisons across offsets (matching setting vs --pf 1)
    fb = res['far_brain']
    for i, nm in ((0, 'full Fourier: 6/8 vs 1'), (1, '7/8 data: 7/8 vs 1'), (2, '6/8 data: 6/8 vs 1')):
        col = {0: 3, 1: 2, 2: 3}[i]
        ratio = fb[:, i, col] / fb[:, i, 1]
        print(f'far_brain {nm}: ratio mean {ratio.mean():.3f} (min {ratio.min():.3f}, max {ratio.max():.3f})')
    if a.json:
        json.dump({m: v.tolist() for m, v in summary.items()}, open(a.json, 'w'))
    if a.figures:
        figures(a, Setup(a.n, a.r, 0, a.phase), a.exe, a.out, summary)


if __name__ == '__main__':
    main()
