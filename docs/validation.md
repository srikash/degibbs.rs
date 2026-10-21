# Validation on a Gibbs phantom

The script is [`validation/phantom.py`](../validation/phantom.py). It regenerates every table and figure below:

```
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --figures docs   # tables and figures
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --phase none     # zero phase
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --phase fast     # phase turns over in 10 voxels
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --side lo        # other end of k-space unsampled
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --snr 20 --offsets 4 --seeds 3
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --snr 8 --offsets 4 --seeds 3
```

#### Setup

- The phantom is the modified Shepp-Logan head (Toft, 1996), scaled so the skull is 100 and the brain 20, with a slowly varying phase map as in real data.
- A 96 × 96 slice is sampled from the object at 8× the resolution, so k-space truncation produces true Gibbs ringing.
- For partial Fourier, the missing lines at one end of k-space are zeroed.
- The magnitude is then taken, as for zero-filled scanner data.
- The truth is the high-resolution object averaged into voxels.
- Errors are RMSE against that truth, averaged over 8 sub-voxel offsets of the object.
- Far-field error is measured more than 4 voxels from any edge, so it shows ringing rather than blur. It is reported for brain tissue (550 pixels, inside the skull with ventricles excluded) and for the background outside the head (3096 pixels).
- Near-edge error is measured within 2 voxels of an edge inside the brain (1952 pixels). The skull ring and its 2-voxel margin (1532 pixels) are reported separately.
- Noise is quoted as brain intensity divided by the standard deviation of each complex channel. It is added on the sampled k-space lines only, so zero-filled partial-Fourier images carry less noise than full-Fourier ones.

#### Results

Far-field RMSE in the brain (lower is better; bold marks the best correction in each row):

| acquired | none | `--pf 1` | `--pf 7/8` | `--pf 6/8` | `--pf 5/8` |
|---|---|---|---|---|---|
| full Fourier | 1.266 | **0.182** | 0.284 | 0.187 | 0.657 |
| 7/8 | 1.229 | 0.559 | **0.251** | 0.576 | 0.714 |
| 6/8 | 1.354 | 0.875 | 0.783 | **0.471** | 0.976 |
| 5/8 | 1.840 | 1.631 | 1.651 | 1.557 | **0.780** |

Far-field RMSE in the background (truth is zero):

| acquired | none | `--pf 1` | `--pf 7/8` | `--pf 6/8` | `--pf 5/8` |
|---|---|---|---|---|---|
| full Fourier | 0.731 | 0.717 | 0.712 | **0.706** | 0.955 |
| 7/8 | 0.753 | 0.740 | 0.737 | **0.731** | 0.936 |
| 6/8 | 0.812 | 0.801 | 0.798 | **0.790** | 0.946 |
| 5/8 | 1.338 | 1.310 | 1.314 | **1.306** | 1.360 |

#### Findings

- In the brain, the matching setting removes 86% (full Fourier), 80% (7/8), 65% (6/8) and 58% (5/8) of the far-field error. It is the best setting in every row with this slowly varying phase.
- On 7/8 and 6/8 data the matching setting has 55% and 46% less error than `--pf 1`. This holds at every offset, between 50% and 61% less for 7/8 and between 45% and 47% less for 6/8.
- On full-Fourier data `--pf 1` beats `--pf 6/8` by 2% to 4% at every offset.
- The `5/8` setting is the worst correction on every other data set, but it still beats no correction inside the brain.
- Near edges inside the brain, the matching setting also lowers the error: 1.37 to 0.92 (full Fourier), 1.27 to 1.05 (7/8), 1.49 to 1.24 (6/8) and 2.24 to 2.03 (5/8).
- In the background, no setting removes more than 4% of the error, because the magnitude of the ringing there is not an oscillation about a plateau.
- Around the skull ring, a thin structure at five times the brain intensity, the partial-Fourier corrections raise the error: 6.06 to 6.95 (7/8), 6.46 to 7.53 (6/8) and 9.98 to 13.18 (5/8). Full-Fourier data improves there, from 6.16 to 5.25. Other thin, bright structures may behave the same way. The whole-image RMSE rises accordingly, for example from 2.76 to 3.01 for 7/8 data.
- Zero phase, or zeroing the other end of k-space, gives the same ordering and gains within 5 percentage points.
- With a fast phase map that turns over in about 10 voxels, the uncorrected error rises sharply for 6/8 and 5/8 data (1.77 and 4.00, against 1.35 and 1.84). The matching setting removes 33% (full Fourier), 36% (7/8), 23% (6/8) and 17% (5/8), and it is no longer always the best: `--pf 7/8` beats `--pf 1` on full-Fourier data (0.62 against 0.69), `--pf 7/8` edges out `--pf 6/8` on 6/8 data (1.34 against 1.36), and `--pf 6/8` beats `--pf 5/8` on 5/8 data (3.23 against 3.34).
- At a noise level of 20, the matching setting removes 51% (full Fourier), 51% (7/8), 48% (6/8) and 52% (5/8) of the error. At 8 the gains are 29%, 34%, 32% and 41%. The truth has no noise, so part of this drop is noise smoothing rather than ringing removal. For example, `--pf 6/8` beats `--pf 1` on full-Fourier data at both noise levels.
- The RPG paper reports robust removal for 7/8 and 6/8 on its own 90 × 90 Shepp-Logan phantom with added phase, and blurring with leftover ringing for 5/8. It attributes the 5/8 result to phase varying within the long ringing interval. The slow phase used here does not test that. The fast phase map does, and there the 5/8 gain falls to 17%, in line with the paper's warning.

#### Figures

![Phantom images and error maps for each acquired partial-Fourier factor](qc_images.png)

*One offset, no noise. Left: truth, uncorrected, corrected with the matching factor and corrected as full Fourier. Right: absolute error maps for the same three images, with the far-field RMSE in the brain and in the background in each title.*

![Phase-encode profiles and RMSE heatmaps](qc_profiles.png)

*Phase-encode profiles through the brain (column 58, one offset), and far-field RMSE heatmaps for the brain and the background (average over 8 offsets).*

This is one phantom. Real data will differ, and the ordering depends on how fast the phase varies. Read the tables for the ordering, not for the size of the gains.

## Complex input

Complex images are treated as full-Fourier data (POCS-reconstructed), so this compares the two paths on the full-Fourier phantom: the magnitude of the image through the magnitude path (`--pf 1`), and the complex image through the complex path. Errors are RMSE of the output magnitude against the truth, averaged over 8 sub-voxel offsets.

| region | none | magnitude path | complex path |
|---|---|---|---|
| far field, brain | 1.266 | 0.182 | 0.185 |
| far field, background | 0.731 | 0.717 | 0.199 |
| near edges, brain | 1.373 | 0.920 | 0.692 |
| skull ring | 6.156 | 5.254 | 3.901 |
| whole image | 2.806 | 2.335 | 1.639 |

- Inside the brain the two paths are equal. The complex path wins elsewhere. In the background the magnitude path removes only 2% of the error, because taking the magnitude turns the ringing into a rectified signal; the complex path removes 73%.
- With noise, the gain shrinks but keeps its direction. At a brain-to-noise ratio of 20 the whole-image error is 2.61 (magnitude path) against 1.92 (complex path); at 8 it is 3.72 against 2.97. Noise is quoted per complex channel, and the truth has no noise, so part of the gain is smoothing.
- The shift for each sample is chosen from the modulus of the complex difference, so the result does not depend on the global phase. A rotation of the whole image by 30, 45 or 90 degrees changes the output by about 1e-8 (RMS, relative to the peak). The alternative of summing the absolute differences of the real and imaginary parts, which the DESIGNER-v2 module uses, changed the output by about 0.3% of the peak at 30 and 45 degrees. Its far-field brain error varied from 0.175 to 0.189 with the rotation, while the modulus gave 0.179 at every rotation. The case for the modulus is that the result does not depend on the arbitrary overall phase, not that it is more accurate.
- Against the DESIGNER-v2 module at `--pf 1`, a modified build that used its measure agreed to 9e-7 (relative RMS), checked once. The shipped measure differs from the module by about 1%, by design.

```
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --complex
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --complex --snr 20 --offsets 4 --seeds 3
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --complex --snr 8 --offsets 4 --seeds 3
```

## POCS and complex input for partial-Fourier data

This answers the question the complex option exists for: acquire partial Fourier, then either zero-fill and use the matching `--pf` on the magnitude, or reconstruct with POCS and use the complex path. The POCS here is basic: the phase comes from the symmetric central lines, followed by 15 iterations. Errors are RMSE against the truth, averaged over 4 sub-voxel offsets, with the slowly varying phase used elsewhere in this document.

Routes: A is the zero-filled magnitude with the matching `--pf` (the usual advice). B is the zero-filled complex image through the complex path, which is the mistake the help text warns about. C is the POCS complex image through the complex path, which is the intended use. D is the POCS magnitude with `--pf 1`.

| acquired | region | none | A | B | C | D |
|---|---|---|---|---|---|---|
| 7/8 | brain, far field | 1.319 | 0.259 | 0.563 | 0.280 | **0.257** |
| 7/8 | background, far field | 0.693 | 0.678 | 0.288 | **0.224** | 0.642 |
| 7/8 | whole image | 2.767 | 2.984 | 2.303 | **1.707** | 2.269 |
| 6/8 | brain, far field | 1.430 | 0.466 | 0.891 | **0.379** | 0.421 |
| 6/8 | background, far field | 0.770 | 0.750 | 0.438 | **0.259** | 0.643 |
| 6/8 | whole image | 2.970 | 3.259 | 3.028 | **1.791** | 2.247 |
| 5/8 | brain, far field | 1.893 | 0.788 | 1.651 | **0.629** | 0.630 |
| 5/8 | background, far field | 1.314 | 1.341 | 1.099 | **0.299** | 0.627 |
| 5/8 | whole image | 4.571 | 5.791 | 5.144 | **1.946** | 2.344 |

- With a slowly varying phase, POCS followed by the complex path is best or tied almost everywhere, and it avoids the skull-ring and whole-image penalty of the partial-Fourier settings. On 5/8 data it cuts the whole-image error from 5.79 (route A) to 1.95.
- Route B gives up most of the benefit. Brain error is about twice that of route A on 7/8 and 6/8 data, and the whole-image error is no better than uncorrected on 6/8 and 5/8 data. This is why complex input that was never POCS-reconstructed should be given as magnitude.
- The comparison depends on the quality of the phase estimate. With the fast phase map (turning over in about 10 voxels), POCS has no good phase to work from. Brain error for route C is then worse than for route A (7/8: 1.218 against 0.727; 6/8: 1.818 against 1.332; 5/8: 3.761 against 3.305) and for 7/8 and 6/8 data it is above the uncorrected value (1.112 and 1.745). Route C still removes most of the background error (0.314 against 0.947 at 7/8).
- This is one phantom and one simple POCS. A scanner's own POCS or homodyne may behave differently.

```
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --pocs --offsets 4
python3 validation/phantom.py --exe target/release/degibbs --out /tmp/pf --pocs --offsets 4 --phase fast
```

## Agreement with the DESIGNER-v2 RPG module

`validation/compare_reference.py` runs this tool and the compiled DESIGNER-v2 RPG Python module on the phantom, for every acquired and correction factor. For `--pf 1`, `7/8` and `6/8` the relative RMS difference is 2e-5 or smaller. For `5/8` magnitude input that module discards the sub-image step and only unrings along the phase-encode axis, while the NYU MATLAB code and this tool run the pipeline described in the paper, so the two differ there.

```
python3 validation/compare_reference.py --exe target/release/degibbs --rpg-path /dir/with/rpg.so --out /tmp/cmp
```

[Back to the README](../README.md)
