"""FSL operations used to build the parcel time series (FSL 6.0.7).

Step 1 (from the images):
  parcel_means   applywarp (trilinear, FNIRT warp of FEAT) to MNI152 2 mm,
                 then fslmeants --label on the Schaefer atlas
  mask_means     fslmeants -m (native-space masks)

Step 2 (from the parcel means):
  glm_resid      fsl_glm --demean --out_res
  bptf_hp        fslmaths -bptf <sigma> -1 (high-pass only)
  denoise        regression, then high-pass, then the temporal mean added
                 back with fslmaths -Tmean / -add

In step 2 the parcel-mean time series are passed to FSL as a 4D image with
one voxel per parcel. Regression and filtering act on each voxel's time
series separately and are linear, so this gives the same result as applying
them to every voxel and averaging afterwards. step2_denoise.py checks this
against FSL's own voxel-wise processing of the ICA-AROMA output
(fmri_preproc).

FSLDIR must point to the FSL installation.
"""
import os
import subprocess
import tempfile

import nibabel as nib
import numpy as np

FSLDIR = os.environ.get('FSLDIR', '')
if not FSLDIR or not os.path.isdir(FSLDIR):
    raise RuntimeError('FSLDIR is not set to an FSL installation (FSL 6.0.7 was used; see README.md).')
_BIN = os.path.join(FSLDIR, 'share', 'fsl', 'bin')
if not os.path.isfile(os.path.join(_BIN, 'fslmaths')):
    _BIN = os.path.join(FSLDIR, 'bin')
_ENV = dict(os.environ, FSLOUTPUTTYPE='NIFTI')


def version():
    try:
        with open(os.path.join(FSLDIR, 'etc', 'fslversion')) as f:
            return f.read().split(':')[0].strip()
    except OSError:
        return 'unknown'


def fsl(tool, *args):
    """Run an FSL tool; stop with its error message if it fails."""
    r = subprocess.run([os.path.join(_BIN, tool), *map(str, args)], env=_ENV,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'{tool} failed ({r.returncode}):\n{r.stdout}\n{r.stderr}')
    return r.stdout


# ---- step 1 -----------------------------------------------------------------

def parcel_means(img, ref, warp, atlas, workdir, interp='trilinear'):
    """Parcel means (volumes x parcels) of a native-space image after
    resampling to standard space: applywarp, then fslmeants --label. Works
    for 3D images too (one row)."""
    std = os.path.join(workdir, 'std.nii')
    txt = os.path.join(workdir, 'meants.txt')
    fsl('applywarp', '-i', img, '-r', ref, '-w', warp, '-o', std, f'--interp={interp}')
    fsl('fslmeants', '-i', std, f'--label={atlas}', '-o', txt)
    os.remove(std)
    return np.atleast_2d(np.loadtxt(txt))


def mask_means(img, mask, workdir):
    """Mean time series of a native-space image within a mask: fslmeants -m."""
    txt = os.path.join(workdir, 'meants.txt')
    fsl('fslmeants', '-i', img, '-m', mask, '-o', txt)
    return np.loadtxt(txt)


# ---- step 2 -----------------------------------------------------------------

def _save(x, path, tr):
    x = np.asarray(x, float)
    T, n = x.shape
    img = nib.Nifti1Image(x.T.reshape(n, 1, 1, T), np.diag([-1.0, 1.0, 1.0, 1.0]))
    img.header.set_data_dtype(np.float64)
    img.header.set_xyzt_units('mm', 'sec')
    img.header['pixdim'][4] = tr
    nib.save(img, path)


def _load(path):
    a = np.asarray(nib.load(path).dataobj, dtype=float)
    return a.reshape(-1, a.shape[-1]).T


class Series:
    """Temporary folder in which time series (volumes x columns) are passed
    to FSL as 4D images with one voxel per column."""

    def __init__(self, tr):
        self.tr = tr
        self._tmp = tempfile.TemporaryDirectory(prefix='fsl_series_')
        self.dir = self._tmp.name
        self._k = 0

    def close(self):
        self._tmp.cleanup()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _path(self, stem):
        self._k += 1
        return os.path.join(self.dir, f'{stem}{self._k}.nii')

    def put(self, x):
        p = self._path('in')
        _save(x, p, self.tr)
        return p

    def glm_resid(self, x_img, Z):
        """fsl_glm -d Z --demean --out_res: residuals of the regression of
        the (demeaned) data on the (demeaned) design."""
        d = os.path.join(self.dir, f'design{self._k}.txt')
        np.savetxt(d, np.asarray(Z, float), fmt='%.17g')
        out = self._path('res')
        fsl('fsl_glm', '-i', x_img, '-d', d, '--demean', f'--out_res={out}')
        return out

    def tmean(self, x_img):
        out = self._path('mean')
        fsl('fslmaths', x_img, '-Tmean', out)
        return out

    def bptf_hp(self, x_img, sigma, add=None):
        """fslmaths -bptf <sigma> -1 (high-pass only; the output has zero
        mean), optionally followed by -add <add>. Computed in double
        precision (-dt double, -odt double)."""
        out = self._path('hp')
        extra = ['-add', add] if add is not None else []
        fsl('fslmaths', '-dt', 'double', x_img, '-bptf', f'{sigma:.10g}', '-1', *extra,
            out, '-odt', 'double')
        return out

    def demean(self, x_img):
        m = self.tmean(x_img)
        out = self._path('dm')
        fsl('fslmaths', '-dt', 'double', x_img, '-sub', m, out, '-odt', 'double')
        return out, m

    def get(self, img):
        return _load(img)

    def denoise(self, x, Z=None, sigma=None):
        """Nuisance regression (if Z is given), then high-pass filter (if
        sigma is given), then the temporal mean of the input added back.
        The series are demeaned (fslmaths, double precision) before
        fsl_glm, because fsl_glm reads its input in single precision."""
        y, m = self.demean(self.put(x))
        if Z is not None:
            y = self.glm_resid(y, Z)
        if sigma is not None:
            return self.get(self.bptf_hp(y, sigma, add=m))
        return self.add(y, m)

    def add(self, y_img, m_img):
        """fslmaths y -add m; returns the series."""
        out = self._path('add')
        fsl('fslmaths', '-dt', 'double', y_img, '-add', m_img, out, '-odt', 'double')
        return self.get(out)


def friston24(par):
    """Friston et al. (1996) 24-parameter model from the six realignment
    parameters p_t: [p_t, p_t^2, p_(t-1), p_(t-1)^2]."""
    p = np.asarray(par, float)[:, :6]
    lag = np.vstack([np.zeros((1, 6)), p[:-1]])
    return np.hstack([p, p ** 2, lag, lag ** 2])
