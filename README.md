# Plasma L-kynurenine and resting-state fMRI Hurst exponents in episodic migraine

Analysis code and coded data for:

> Szuly A, Tuka B, Szabó N, Faragó P, Csomós M, Kocsis K, Körtési T, Vécsei L, Pereira JB, Veréb D, Kincses ZT. Plasma L-kynurenine levels are associated with resting-state cortical dynamics in episodic migraine.

The scripts in this repository produce the numbers, figures (Figures 1–3, S1) and supplementary tables (Tables S1–S8) of the article. The random seed is 20261105.

Archived on Zenodo: version 1.0.0, used for the article, https://doi.org/10.5281/zenodo.24278680; all versions, https://doi.org/10.5281/zenodo.24278679.

## Contents

| Path | Content |
|---|---|
| `run_all.sh` | runs the whole analysis |
| `step0_export_inputs.m` … `step5_figures.R` | the six steps of the pipeline (see below) |
| `check_hurst.m` | recomputes the per-participant Hurst exponents from the deposited time series and compares them with `results/`; needs no clinical data |
| `fsl_tools.py` | the FSL calls of steps 1 and 2 (`applywarp`, `fslmeants`, `fsl_glm`, `fslmaths`) |
| `figure1_pathway.py` | Figure 1 (schematic of the kynurenine pathway) |
| `lib/` | MATLAB functions: Hurst estimator, partial Spearman correlation, permutation, bootstrap, Benjamini–Hochberg correction |
| `data/` | coded parcel time series and head-motion summaries of the 21 participants |
| `results/` | all results in long format, per-participant Hurst exponents, run logs |
| `tables/` | Supplementary Tables S1–S8 (comma-separated, unrounded values) |
| `figures/` | Figures 1–3 and S1 (PDF and 600 dpi PNG) |
| `toolboxes/README.md` | how to obtain the third-party MATLAB toolboxes |

## Quick start

```bash
./run_all.sh
```

This runs steps 2 to 5 from the files in `data/`. Steps 3 to 5 need the clinical table `data/participants.csv`, which cannot be shared (see [Clinical table](#clinical-table)). Without it, the script runs step 2 and `check_hurst.m`, and stops there. `check_hurst.m` recomputes the Hurst exponents of all participants and parcels from the output of step 2, compares them with `results/hurst_parcels.csv` and `results/hurst_networks.csv` and prints the largest difference; with the deposited files this is 5 × 10⁻¹⁶, which comes from the rounding of the values in the CSV files. It stops with an error if any value differs by more than 10⁻¹². It writes no files and can also be run separately (`matlab -batch check_hurst`). Step 2 calls FSL, which is found through the environment variable `FSLDIR` (set by the standard FSL installation). The interpreters can be set with the environment variables `PYTHON`, `MATLAB` and `RSCRIPT`, and the folder of the MATLAB toolboxes with `HURST_TOOLBOXES` (default: `toolboxes/`).

## Requirements

- MATLAB (R2024a was used) with the Statistics and Machine Learning Toolbox, and the nonfractal and WMTSA toolboxes; see `toolboxes/README.md`. The toolboxes are not redistributed here.
- FSL (steps 1 and 2; FSL 6.0.7 was used): `applywarp` and `fslmeants` in step 1, `fsl_glm` and `fslmaths` in step 2. `FSLDIR` must point to the installation.
- Python 3 with numpy, scipy and nibabel (steps 1 and 2; Python 3.13, numpy 2.3, scipy 1.16, nibabel 5.4 were used) and matplotlib (Figure 1).
- R with ggplot2, dplyr, readr, patchwork and scales (step 5; R 4.6, ggplot2 4.0, dplyr 1.2, readr 2.2, patchwork 1.3, scales 1.4 were used).

## Pipeline

| Step | Script | Input | Output |
|---|---|---|---|
| 0 | `step0_export_inputs.m` | original analysis workspace and official Schaefer label file | `data/participants.csv`, `data/parcels.csv`, `data/original_submission_ts.mat`, `private/subject_key.csv` |
| 1 | `step1_extract.py` | image folders, `private/subject_key.csv`, Schaefer atlas | `data/extracted/S##.npz` |
| 2 | `step2_denoise.py` | `data/extracted/` | `data/timeseries/S##.csv` (primary series), `data/denoised_variants.mat`, `data/motion.csv`, `data/parcel_coverage.csv` |
| – | `check_hurst.m` | `data/denoised_variants.mat`, `data/parcels.csv`, `results/hurst_parcels.csv`, `results/hurst_networks.csv` | none (prints the comparison) |
| 3 | `step3_main_analysis.m` | `data/` | `results/*.csv`, Tables S1, S2, S6, S7 |
| 4 | `step4_sensitivity.m` | `data/`, `results/step3_workspace.mat` (written by step 3) | `results/*.csv`, Tables S3, S4, S5, S8 |
| 5 | `step5_figures.R` | `results/` | Figures 2, 3 and S1, `figures/captions.md` |

Steps 0 and 1 need the images and a key linking the participant codes to the image folders, which are not shared. They are included to document how the files in `data/` were made, and their paths refer to the folder layout of the original project. `./run_all.sh --from-images` runs steps 0 to 5 where these files are available. Figure 1 is drawn by `python figure1_pathway.py`.

## Data

Participants are identified by the codes S01–S21.

| File | Content |
|---|---|
| `data/extracted/S##.npz` | Output of step 1 for each participant (198 volumes, 100 parcels): `feat`, parcel means of the FEAT pre-statistics output; `aroma`, parcel means of the ICA-AROMA (non-aggressive) output; `fsl_preproc`, parcel means of FSL's own fully pre-processed series (used to check step 2, and as a sensitivity variant); `motion`, the six MCFLIRT realignment parameters (rotations in rad, translations in mm); `fd_rel`, MCFLIRT relative root-mean-square displacement (mm); `wmcsf_feat` and `wmcsf_aroma`, mean white-matter and CSF signals; `mix_motion`, time courses of the components that ICA-AROMA classified as motion; `n_ics`, number of ICA components; `coverage`, share of each parcel inside the brain mask |
| `data/timeseries/S##.csv` | Primary parcel time series, 198 volumes × 100 parcels (columns `P001`–`P100`): Friston 24-parameter regression, high-pass filter at 0.01 Hz, temporal mean added back |
| `data/denoised_variants.mat` | All ten denoising variants of Supplementary Table S4, each 198 × 100 × 21 (`f24_hp` is the primary series), with the variant names (`variants`) and participant codes (`ids`) |
| `data/original_submission_ts.mat` | `ts`, 198 × 100 × 21: the parcel time series analysed in the first submitted version of the article, which were extracted before this pipeline was set up (last row of Table S4) |
| `data/parcels.csv` | Schaefer-100 parcel labels and 7-network assignment from the official label file; `network_original_submission` is the assignment used in the first submitted version, which differed for parcels 79–83, 88 and 89 |
| `data/motion.csv` | Per participant: mean and maximum framewise displacement (Jenkinson and Power definitions), percentage of volumes with FD (Power) > 0.5 mm, maximum translation and rotation, number of ICA-AROMA motion components and of all components |
| `data/parcel_coverage.csv` | Share of each parcel inside the brain mask, per participant |

Variants in `data/denoised_variants.mat`, in the row order of Table S4: `f24_hp`, `hp_then_f24`, `f24deriv_hp`, `f24_wmcsf_hp`, `mp6_hp`, `mp6_wmcsf_hp`, `mcf_hp` (high-pass filter only), `aroma_nonaggr_wmcsf_hp`, `aroma_aggr_hp`, `original_submission`.

### Clinical table

`data/participants.csv` holds the age, aura status and plasma metabolite concentrations of the participants. It cannot be shared, because the ethics approval of the study does not cover the sharing of individual clinical data. Its format is given in `data/participants_template.csv`. With the clinical table in `data/`, `./run_all.sh` reproduces all results. From the deposited data alone, step 2 (the denoising of the parcel time series) and the estimation of the Hurst exponents (`check_hurst.m`) can be repeated; the per-participant Hurst exponents and all results are provided in `results/` and `tables/`.

Columns of the clinical table:

| Column | Content |
|---|---|
| `id` | participant code, S01–S21 |
| `aura` | 1 = migraine with aura, 0 = migraine without aura |
| `age` | years |
| `TRP`, `L_KYN`, `KYNA`, `ANA`, `HK3`, `XA`, `HANA3`, `PICA`, `QUINA` | plasma concentrations (nmol/L) of tryptophan, L-kynurenine, kynurenic acid, anthranilic acid, 3-hydroxykynurenine, xanthurenic acid, 3-hydroxyanthranilic acid, picolinic acid and quinolinic acid |

Step 3 also writes `results/fig_subjects.csv` (L-KYN, aura status and SMN Hurst exponent per participant, for the figures); it is excluded from the repository for the same reason (see `.gitignore`).

## Preprocessing

The starting point is the FEAT pre-statistics output of each participant: first two volumes removed, motion correction with MCFLIRT, brain extraction with BET, spatial smoothing with a 5 mm FWHM kernel, grand-mean scaling, no temporal filtering; 198 volumes remain (TR = 3 s).

Primary pipeline (`f24_hp`), as described in the Methods of the article:

1. Parcel means: the native-space data are resampled to MNI152 2 mm space with the FNIRT warp of FEAT (`applywarp`, trilinear interpolation) and averaged within the 100 parcels of the Schaefer atlas, 7-network version (`fslmeants --label`).
2. Regression of the Friston 24-parameter motion model: 6 realignment parameters, their values at the previous volume, and the squares of both (`fsl_glm --demean --out_res`).
3. High-pass filter at 0.01 Hz (`fslmaths -bptf 16.67 -1`; sigma = 1/(2 × 0.01 Hz × TR) = 16.67 volumes).
4. The temporal mean is added back (`fslmaths -Tmean`, `-add`).

All four are FSL 6.0.7 operations, called from `fsl_tools.py`. Steps 2 to 4 act on each voxel's time series separately and are linear, so applying them to parcel means is equivalent to applying them to every voxel and averaging afterwards. In step 2 the parcel-mean time series are therefore passed to FSL as a 4D image with one voxel per parcel. They are demeaned in double precision before `fsl_glm`, which reads its input in single precision. The other variants in step 2 are sensitivity analyses (Table S4).

Step 2 checks this equivalence on every run with the ICA-AROMA output. For this output, the regression of the white-matter and CSF signals and the high-pass filter were applied voxel-wise by FSL in the original analysis (`fmri_preproc`). Applying the same operations to the parcel means gives the parcel means of the voxel-wise result with a relative root-mean-square error of 2.8 × 10⁻⁶ (median over 21 participants × 100 parcels; worst parcel 1.2 × 10⁻⁵).

## Hurst exponent

Main analysis: `bfn_mfin_ml` of the nonfractal toolbox with its default options: discrete wavelet transform with the LA8 filter and periodic boundary handling, boundary coefficients removed, the first 128 of the 198 volumes (the largest power of two; no padding), wavelet levels 1 to 4 (0.010–0.167 Hz), maximum-likelihood estimation of d = H − 0.5 within [−0.5, 0.5].

`lib/hurst_wavelet_ml.m` is an independent implementation of the same estimator. With the options `'dwt', 'la8', 'variance'` it gives the same values as the toolbox (maximum absolute difference 0; checked in step 3 and in `check_hurst.m`). It also provides the alternatives of Table S3: a maximal-overlap discrete wavelet transform (`'modwt'`), which is not restricted to a power-of-two length and can use all 198 volumes, and a profile-likelihood variant.

Network Hurst exponents are the mean over the parcels of a network; the whole-cortex value is the mean over all 100 parcels.

## Statistics

- Partial Spearman correlation between a plasma measure and a Hurst exponent, adjusted for age and aura status (`lib/ps_corr.m`).
- Permutation p-value: the plasma measure is permuted across participants, 20,000 permutations, p = (b + 1)/(m + 1) (`lib/ps_perm.m`).
- 95 % confidence interval: bias-corrected and accelerated bootstrap, 10,000 resamples drawn within aura strata (`lib/ps_boot.m`, `lib/bca.m`).
- Benjamini–Hochberg correction within the primary family (L-KYN × 7 networks), within the secondary family (9 other measures × 7 networks) and across all 70 tests (`lib/bh_fdr.m`).
- All permutation and bootstrap indices are drawn once by `lib/make_resamples.m` from `rng(20261105, 'twister')`. Steps 3 and 4 use the same resamples, so the main-analysis row of each sensitivity table is identical to the main result. Two subset analyses draw their own indices: the participants without aura (inside `make_resamples`) and the exclusion of the two participants with more than 5 % of volumes above 0.5 mm framewise displacement (seed 20261106, end of step 4).

## Tables

| File | Content |
|---|---|
| `tables/Table_S1_network_by_measure.csv` | all 70 tests: rho, 95 % CI, permutation p, q within family, q across 70 |
| `tables/Table_S2_leave_one_out.csv` | leave-one-out estimates, 7 networks and whole cortex |
| `tables/Table_S3_wavelet_settings.csv` | 12 estimator settings |
| `tables/Table_S4_denoising_variants.csv` | 10 denoising variants, including ICA-AROMA and the series of the first submitted version |
| `tables/Table_S5_window_position.csv` | 128-volume windows at 8 positions, two estimators |
| `tables/Table_S6_head_motion.csv` | head motion against L-KYN and Hurst exponents; framewise displacement as a covariate |
| `tables/Table_S7_metabolite_ratios.csv` | four metabolite ratios × 7 networks |
| `tables/Table_S8_additional_checks.csv` | ICA-AROMA motion-component share as covariate; poorly covered parcels excluded; two participants with more than 5 % of volumes above 0.5 mm framewise displacement excluded |

`results/` holds the same results in long format, the per-participant Hurst exponents (`hurst_parcels.csv`, `hurst_networks.csv`) and the logs of steps 3 and 4.

## Licence

Code: MIT licence (`LICENSE`). Data, results, tables and figures (`data/`, `results/`, `tables/`, `figures/`): Creative Commons Attribution 4.0 International (CC BY 4.0; `LICENSE-DATA.md`).

## Contact

Attila Szuly, Department of Radiology, Albert Szent-Györgyi Health Centre, University of Szeged, Szeged, Hungary (szuly.attila.krisztian@szte.hu)
