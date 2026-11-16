"""step1_extract.py

Extracts the data needed by the later steps from the image folder of each
participant and saves them under the participant code (S01-S21). Parcel
means: each native-space image is resampled to MNI152 2 mm with the FNIRT
warp of FEAT (FSL applywarp, trilinear) and averaged within the parcels of
the Schaefer atlas (FSL fslmeants --label); see fsl_tools.py.

  feat         parcel means of the FEAT pre-statistics output
               (rsn.ica/orig_filtered_func_data: first 2 volumes removed,
               MCFLIRT, BET, 5 mm SUSAN smoothing, grand-mean scaling, no
               temporal filtering), 198 x 100
  aroma        parcel means of the ICA-AROMA (non-aggressive) output
  fsl_preproc  parcel means of FSL's own fmri_preproc (AROMA + WM/CSF
               regression + high-pass, voxel-wise), used to check step 2 and as
               a sensitivity variant
  motion       six MCFLIRT realignment parameters (rot x,y,z in rad; trans
               x,y,z in mm)
  fd_rel       MCFLIRT relative RMS displacement per volume pair (mm)
  wmcsf_feat   mean white-matter and CSF signal of the FEAT output (on-disk
               native-space masks, thresholded at 0.9; fslmeants -m)
  wmcsf_aroma  the same from the AROMA output (nuisance_design.txt on disk)
  mix_motion   time courses of the ICA components AROMA classified as motion
  n_ics        number of ICA components
  coverage     share of each parcel that falls inside the FEAT brain mask
               (the mask resampled in the same way, then averaged)

Needs FSL (FSLDIR), the image folders and private/subject_key.csv.
The resampled images are written to a temporary folder and deleted after use.
Usage:

    python -I step1_extract.py [NIFTI_DIR]
"""
import csv
import os
import sys
import tempfile

import nibabel as nib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fsl_tools import mask_means, parcel_means, version  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(HERE))
NIFTI = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'migraine_kyn', 'niftis')
ATLAS = os.path.join(ROOT, 'toolboxes', 'schaefer',
                     'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
OUT = os.path.join(HERE, 'data', 'extracted')
os.makedirs(OUT, exist_ok=True)


at = nib.load(ATLAS)
key = list(csv.DictReader(open(os.path.join(HERE, 'private', 'subject_key.csv'))))
assert len(key) == 21
print(f'FSL {version()}', flush=True)
tmp = tempfile.TemporaryDirectory(prefix='step1_')     # removed on exit, also after an error
work = tmp.name
for row in key:
    sid, d = row['id'], os.path.join(NIFTI, row['folder'])
    ica, aro, nui = os.path.join(d, 'rsn.ica'), os.path.join(d, 'ica_aroma'), os.path.join(d, 'nuisance')
    ref = os.path.join(ica, 'reg', 'standard.nii.gz')
    warp = os.path.join(ica, 'reg', 'example_func2standard_warp.nii.gz')
    r = nib.load(ref)
    assert r.shape == at.shape and np.allclose(r.affine, at.affine), 'reference and atlas grids differ'
    feat_img = os.path.join(ica, 'orig_filtered_func_data.nii.gz')
    feat = parcel_means(feat_img, ref, warp, ATLAS, work)
    aroma = parcel_means(os.path.join(aro, 'denoised_func_data_nonaggr.nii.gz'), ref, warp, ATLAS, work)
    fsl_preproc = parcel_means(os.path.join(d, 'fmri_preproc.nii.gz'), ref, warp, ATLAS, work)
    motion = np.loadtxt(os.path.join(ica, 'mc', 'prefiltered_func_data_mcf.par'))[:, :6]
    fd_rel = np.loadtxt(os.path.join(ica, 'mc', 'prefiltered_func_data_mcf_rel.rms'))
    fd_rel_mean = float(np.loadtxt(os.path.join(ica, 'mc', 'prefiltered_func_data_mcf_rel_mean.rms')))
    wmcsf_feat = np.column_stack([mask_means(feat_img, os.path.join(nui, 'wm_native.nii.gz'), work),
                                  mask_means(feat_img, os.path.join(nui, 'csf_native.nii.gz'), work)])
    wmcsf_aroma = np.loadtxt(os.path.join(nui, 'nuisance_design.txt'))
    mix = np.loadtxt(os.path.join(aro, 'melodic.ica', 'melodic_mix'))
    mot = np.atleast_1d(np.loadtxt(os.path.join(aro, 'classified_motion_ICs.txt'),
                                   delimiter=',')).astype(int) - 1
    coverage = parcel_means(os.path.join(ica, 'mask.nii.gz'), ref, warp, ATLAS, work).ravel()

    T = feat.shape[0]
    assert T == 198 and feat.shape == aroma.shape == fsl_preproc.shape == (198, 100) and coverage.shape == (100,)
    assert motion.shape == (T, 6) and wmcsf_feat.shape == wmcsf_aroma.shape == (T, 2)
    assert fd_rel.shape == (T - 1,) and mix.shape[0] == T
    assert abs(fd_rel.mean() - fd_rel_mean) < 1e-4
    np.savez_compressed(os.path.join(OUT, f'{sid}.npz'), feat=feat, aroma=aroma,
                        fsl_preproc=fsl_preproc, motion=motion, fd_rel=fd_rel,
                        wmcsf_feat=wmcsf_feat, wmcsf_aroma=wmcsf_aroma,
                        mix_motion=mix[:, mot], n_ics=mix.shape[1], coverage=coverage)
    print(f'{sid}: {T} volumes, {len(mot)}/{mix.shape[1]} motion components, '
          f'mean FD {fd_rel_mean:.3f} mm, parcel coverage min {coverage.min():.2f}, '
          f'parcels < 0.5 covered: {int((coverage < 0.5).sum())}', flush=True)
tmp.cleanup()
print('done')
