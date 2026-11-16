# Third-party toolboxes

Steps 3 and 4 and `check_hurst.m` need two MATLAB toolboxes, which are not redistributed here. `lib/setup_paths.m` looks for them in this folder (or in the folder named by the environment variable `HURST_TOOLBOXES`) with this layout:

```
toolboxes/
  nonfractal/   nonfractal toolbox (contains m/bfn_mfin_ml.m)
  wmtsa/        WMTSA Wavelet Toolkit for MATLAB 0.2.6
  wmtsa_mex/    dwtj and modwtj of WMTSA, compiled with mex
```

| Toolbox | Version used | Source | SHA-256 of the archive |
|---|---|---|---|
| nonfractal | commit `38accd45956af6c8a65a0d0b4005d04261e73b76` | https://github.com/wonsang/nonfractal | `579af96feaed6a7f405ac634b20742a03c6410fec7cb618c93f4f16372e0842a` (GitHub zip of that commit, downloaded 2026-11-05) |
| WMTSA Wavelet Toolkit for MATLAB | 0.2.6 | The original site (atmos.washington.edu/~wmtsa) is offline; Internet Archive copy: `http://web.archive.org/web/20170425151709id_/http://atmos.washington.edu/%7Ewmtsa//wmtsa-matlab-0.2.6.zip` | `15b47ad1e34f7326d748a0dce046c8a58a06faa10c7a68858007292de13ac339` |

## Installation

From this folder:

```bash
git clone https://github.com/wonsang/nonfractal.git
git -C nonfractal checkout 38accd45956af6c8a65a0d0b4005d04261e73b76
curl -L -o wmtsa-matlab-0.2.6.zip "http://web.archive.org/web/20170425151709id_/http://atmos.washington.edu/%7Ewmtsa//wmtsa-matlab-0.2.6.zip"
shasum -a 256 wmtsa-matlab-0.2.6.zip
unzip -q wmtsa-matlab-0.2.6.zip -d wmtsa
mkdir wmtsa_mex
```

WMTSA's `dwt.m` and `modwt.m` call `dwtj` and `modwtj`, which are distributed as C source. Compile them in MATLAB (a C compiler supported by `mex` is needed), from this folder:

```matlab
mex -outdir wmtsa_mex wmtsa/wmtsa-matlab-0.2.6/wmtsa/dwt/dwtj.c
mex -outdir wmtsa_mex wmtsa/wmtsa-matlab-0.2.6/wmtsa/dwt/modwtj.c
```

With these files, `bfn_mfin_ml` with its default options gives the Hurst exponents reported in the article. Step 3 and `check_hurst.m` compare it with `lib/hurst_wavelet_ml.m` and print the maximum absolute difference.

## Schaefer atlas (steps 0 and 1 only)

Schaefer et al. (2018), 100 parcels, 7 networks, FSL MNI152 2 mm, from the CBIG repository (https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI):

- `Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz` (SHA-256 `101cb2e643531826e5d52d35751f3cf7134c394b08b8ce0affa259a224350f67`)
- `freeview_lut/Schaefer2018_100Parcels_7Networks_order.txt` (SHA-256 `ca212790d0103949c8cfa0f7689f5935cad0d1fbfc06506f1fdf7986bee72786`)

The parcel labels and network assignments derived from these files are in `data/parcels.csv`, so steps 2 to 5 do not need the atlas.
