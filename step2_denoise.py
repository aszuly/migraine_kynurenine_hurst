"""step2_denoise.py

Builds the denoised parcel time series from data/extracted/S##.npz. The
regression and filtering are done by FSL (fsl_glm, fslmaths -bptf; see
fsl_tools.py), so FSL is needed (FSLDIR).

Primary pipeline (described in the Methods of the article):

    FEAT pre-statistics output
      -> regression of the Friston 24-parameter motion model
      -> high-pass filter, 0.01 Hz (FSL -bptf, sigma = 1/(2*0.01*TR) volumes)
      -> temporal mean added back

Sensitivity variants (reported in a supplementary table only):

    hp_then_f24      high-pass first (data and regressors), then Friston-24
    f24deriv_hp      24 parameters as 6 + temporal derivatives + squares
    f24_wmcsf_hp     Friston-24 + white-matter and CSF signals
    mp6_hp           six realignment parameters
    mp6_wmcsf_hp     six realignment parameters + WM and CSF signals
    mcf_hp           no nuisance regression (motion correction only)
    aroma_nonaggr_wmcsf_hp   ICA-AROMA (non-aggressive) + WM/CSF + high-pass,
                             as produced by FSL on disk (fmri_preproc,
                             processed voxel-wise)
    aroma_aggr_hp    full regression of the AROMA motion components
    original_submission      series analysed in the original submission

Outputs
    data/denoised_variants.mat   one 198 x 100 x 21 array per variant
    data/timeseries/S##.csv      primary series, volumes x parcels
    data/motion.csv              head-motion summary per participant
    data/parcel_coverage.csv     share of each parcel inside the brain mask

Usage:  python -I step2_denoise.py
"""
import os
import sys

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fsl_tools import Series, friston24, version  # noqa: E402

TR = 3.0
CUTOFF_HZ = 0.01
SIGMA = 1.0 / (2.0 * CUTOFF_HZ * TR)        # 16.667 volumes
IDS = [f'S{k:02d}' for k in range(1, 22)]
PRIMARY = 'f24_hp'
VARIANTS = [PRIMARY, 'hp_then_f24', 'f24deriv_hp', 'f24_wmcsf_hp', 'mp6_hp', 'mp6_wmcsf_hp',
            'mcf_hp', 'aroma_nonaggr_wmcsf_hp', 'aroma_aggr_hp', 'original_submission']


def friston24_deriv(p):
    d = np.vstack([np.zeros((1, 6)), np.diff(p, axis=0)])
    return np.hstack([p, d, p ** 2, d ** 2])


def hp_then_regress(S, x, Z):
    """High-pass filter of the data and of the regressors, then regression,
    then the temporal mean added back."""
    xd, m = S.demean(S.put(x))
    zd, _ = S.demean(S.put(Z))
    r = S.glm_resid(S.bptf_hp(xd, SIGMA), S.get(S.bptf_hp(zd, SIGMA)))
    return S.add(r, m)


print(f'FSL {version()}')
orig = sio.loadmat(os.path.join(HERE, 'data', 'original_submission_ts.mat'))['ts']
P = {v: np.zeros((198, 100, 21)) for v in VARIANTS}
P['original_submission'] = orig.astype(float)
motion_rows = []
val = []
cover = []
for k, sid in enumerate(IDS):
    E = np.load(os.path.join(HERE, 'data', 'extracted', f'{sid}.npz'))
    x, mp = E['feat'].astype(float), E['motion'].astype(float)
    cover.append(E['coverage'])
    wc = E['wmcsf_feat'].astype(float)
    F = friston24(mp)
    with Series(TR) as S:
        P['f24_hp'][:, :, k] = S.denoise(x, F, SIGMA)
        P['hp_then_f24'][:, :, k] = hp_then_regress(S, x, F)
        P['f24deriv_hp'][:, :, k] = S.denoise(x, friston24_deriv(mp), SIGMA)
        P['f24_wmcsf_hp'][:, :, k] = S.denoise(x, np.hstack([F, wc]), SIGMA)
        P['mp6_hp'][:, :, k] = S.denoise(x, mp, SIGMA)
        P['mp6_wmcsf_hp'][:, :, k] = S.denoise(x, np.hstack([mp, wc]), SIGMA)
        P['mcf_hp'][:, :, k] = S.denoise(x, None, SIGMA)
        P['aroma_nonaggr_wmcsf_hp'][:, :, k] = E['fsl_preproc']
        P['aroma_aggr_hp'][:, :, k] = S.denoise(x, E['mix_motion'].astype(float), SIGMA)

        # Check: the same FSL operations applied to the parcel means of the
        # AROMA output must reproduce the parcel means of FSL's voxel-wise
        # processing of that output (fmri_preproc: WM/CSF regression and
        # high-pass filter applied to every voxel).
        rep = S.denoise(E['aroma'].astype(float), E['wmcsf_aroma'].astype(float), SIGMA)
    ref = E['fsl_preproc'].astype(float)
    val.append(np.sqrt(((rep - ref) ** 2).sum(0) / ((ref - ref.mean(0)) ** 2).sum(0)))

    # Head motion
    d = np.abs(np.diff(mp, axis=0))
    fd_power = d[:, 3:6].sum(1) + 50.0 * d[:, 0:3].sum(1)      # rotations on a 50 mm sphere
    fd_rel = E['fd_rel']
    motion_rows.append((sid, fd_rel.mean(), fd_rel.max(), fd_power.mean(), fd_power.max(),
                        100.0 * (fd_power > 0.5).mean(), np.abs(mp[:, 3:6]).max(),
                        np.degrees(np.abs(mp[:, 0:3]).max()),
                        E['mix_motion'].shape[1], int(E['n_ics'])))

    np.savetxt(os.path.join(HERE, 'data', 'timeseries', f'{sid}.csv'), P[PRIMARY][:, :, k],
               delimiter=',', fmt='%.6f', header=','.join(f'P{j:03d}' for j in range(1, 101)),
               comments='')

val = np.array(val)
print(f'parcel-level vs voxel-wise FSL processing (fmri_preproc): relative RMS error '
      f'median {np.median(val):.2e}, worst parcel {val.max():.2e}')
for v in VARIANTS:
    a = P[v]
    print(f'{v:24s} NaN: {bool(np.isnan(a).any())}  zero-variance series: {int((a.std(0) == 0).sum())}')
sio.savemat(os.path.join(HERE, 'data', 'denoised_variants.mat'),
            dict(P, variants=np.array(VARIANTS, dtype=object), primary=PRIMARY,
                 ids=np.array(IDS, dtype=object)), do_compression=True)

with open(os.path.join(HERE, 'data', 'motion.csv'), 'w') as f:
    f.write('id,FD_jenkinson_mean,FD_jenkinson_max,FD_power_mean,FD_power_max,'
            'pct_vols_FDpower_gt_0p5,max_abs_translation_mm,max_abs_rotation_deg,'
            'aroma_motion_ICs,ICs_total\n')
    for r in motion_rows:
        f.write(f'{r[0]},{r[1]:.5f},{r[2]:.5f},{r[3]:.5f},{r[4]:.5f},{r[5]:.3f},{r[6]:.4f},{r[7]:.4f},{r[8]},{r[9]}\n')
with open(os.path.join(HERE, 'data', 'parcel_coverage.csv'), 'w') as f:
    f.write('id,' + ','.join(f'P{j:03d}' for j in range(1, 101)) + '\n')
    for sid, c in zip(IDS, cover):
        f.write(sid + ',' + ','.join(f'{v:.4f}' for v in c) + '\n')
m = np.array([r[1:] for r in motion_rows], float)
print(f'mean FD (Jenkinson): median {np.median(m[:, 0]):.3f} mm, range {m[:, 0].min():.3f}-{m[:, 0].max():.3f}; '
      f'participants > 0.2 mm: {int((m[:, 0] > 0.2).sum())}')
print(f'mean FD (Power):     median {np.median(m[:, 2]):.3f} mm, range {m[:, 2].min():.3f}-{m[:, 2].max():.3f}')
print(f'max absolute translation {m[:, 5].max():.2f} mm, rotation {m[:, 6].max():.2f} deg')
print('done')
