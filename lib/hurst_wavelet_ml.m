function [H, info] = hurst_wavelet_ml(X, transform, filt, sigma, maxlevel)
% HURST_WAVELET_ML  Hurst exponents by wavelet maximum likelihood.
%
%   H = hurst_wavelet_ml(X)
%   H = hurst_wavelet_ml(X, transform, filt, sigma, maxlevel)
%   [H, info] = hurst_wavelet_ml(...)
%
% Fits a fractionally integrated noise (d = H - 0.5, bounded to
% [-0.5, 0.5]) to each column of X by wavelet maximum likelihood, using
% the spectral model of the nonfractal toolbox (bfn_finsdf). For each
% wavelet level j the likelihood uses the sum of squared non-boundary
% coefficients S_j and the number of coefficients N_j:
%       L(d, s2) = sum_j N_j log s_j(d) + sum_j S_j / (s_j(d) s2)
% where s_j(d) is the band-integrated spectral density.
%
% Inputs
%   X          T x Q matrix, one time series per column.
%   transform  'dwt'     The same as bfn_mfin_ml with default options:
%                        series truncated to the largest dyadic length,
%                        DWT with periodic boundary, boundary coefficients
%                        removed, N_j = N / 2^j (the nominal count,
%                        including the removed coefficients).
%              'dwt_eff' As 'dwt', but N_j = number of coefficients
%                        actually summed into S_j.
%              'modwt'   The whole series (any length T). MODWT with
%                        circular boundary; the first L_j - 1 coefficients
%                        of level j, L_j = (2^j - 1)(L - 1) + 1, are
%                        boundary-affected and removed. S_j is the sum of
%                        squares of the M_j = T - L_j + 1 remaining
%                        coefficients and N_j = M_j / 2^j, the equivalent
%                        number of DWT coefficients (Percival & Walden
%                        2000, eq. 314c). Averages the DWT statistics over
%                        all circular shifts of the series.
%              Default 'modwt'.
%   filt       Wavelet filter name, e.g. 'la8' (default) or 'haar'.
%   sigma      How the innovation variance s2 is obtained.
%              'variance' As bfn_mfin_ml with default options
%                         (omegamode 'cov'): s2 is computed from the
%                         sample variance of the series and the current
%                         d, and d and s2 are updated in turn (fin_ml).
%              'profile'  s2 is set to its maximum-likelihood value given
%                         d, s2 = sum_j(S_j / s_j(d)) / sum_j N_j, and the
%                         resulting profile likelihood is minimised over d.
%              Default 'profile'.
%   maxlevel   Highest level fitted (default: all levels). Levels are the
%              WMTSA 'conservative' choice, J0 = floor(log2(N/(L-1) - 1));
%              level j spans 1/2^(j+1) to 1/2^j cycles per sample.
%
% Requires the nonfractal and WMTSA toolboxes on the path, with WMTSA's
% dwtj/modwtj compiled as MEX files (see ../toolboxes/README.md).

if nargin < 2 || isempty(transform), transform = 'modwt'; end
if nargin < 3 || isempty(filt), filt = 'la8'; end
if nargin < 4 || isempty(sigma), sigma = 'profile'; end
if nargin < 5, maxlevel = []; end

[T, Q] = size(X);
switch transform
    case {'dwt', 'dwt_eff'}
        N = 2^floor(log2(T));
        X = X(1:N, :);
    case 'modwt'
        N = T;
    otherwise
        error('hurst_wavelet_ml:transform', 'Unknown transform ''%s''.', transform);
end
wtf = dwt_filter(filt);
L = wtf.L;

H = nan(1, Q);
for q = 1:Q
    x = X(:, q);
    switch transform
        case {'dwt', 'dwt_eff'}
            W = dwt(x, filt, 'conservative', 'periodic');
            W = bfn_dwt_brick_wall(W, filt, N);
            J = numel(W);
            S = zeros(J, 1);
            Nj = zeros(J, 1);
            for j = 1:J
                w = W{j}(~isnan(W{j}));
                S(j) = sum(w .* w);
                if strcmp(transform, 'dwt')
                    Nj(j) = numel(W{j});
                else
                    Nj(j) = numel(w);
                end
            end
        case 'modwt'
            W = modwt(x, filt, 'conservative', 'circular');
            J = size(W, 2);
            S = zeros(J, 1);
            Nj = zeros(J, 1);
            for j = 1:J
                Lj = (2^j - 1) * (L - 1) + 1;
                w = W(Lj:N, j);
                S(j) = sum(w .* w);
                Nj(j) = numel(w) / 2^j;
            end
    end
    if ~isempty(maxlevel)
        assert(maxlevel <= J, 'hurst_wavelet_ml: only %d levels available.', J);
        J = maxlevel;
        S = S(1:J);
        Nj = Nj(1:J);
    end
    switch sigma
        case 'variance'
            H(q) = fin_ml(S, var(x), Nj, J) + 0.5;
        case 'profile'
            H(q) = profile_ml(S, Nj, J) + 0.5;
        otherwise
            error('hurst_wavelet_ml:sigma', 'Unknown sigma mode ''%s''.', sigma);
    end
end

info = struct('transform', transform, 'filter', filt, 'sigma', sigma, ...
              'N', N, 'J', J, ...
              'Nj', Nj', 'band_cycles_per_sample', ...
              [1 ./ 2.^((1:J)' + 1), 1 ./ 2.^(1:J)']);
end

% -------------------------------------------------------------------------
function d = fin_ml(S, va, Nj, J)
% Iterative ML of nonfractal's fin_ml with its default settings
% (omegamode 'cov', init d = 0.1, bounds [-0.5, 0.5], abstol 1e-10,
% maxit 100). Same operations in the same order, so that the 'dwt'
% transform reproduces bfn_mfin_ml.
    Js = 1:J;
    lb = -0.5;
    ub = 0.5;
    abstol = 1e-10;
    maxit = 100;
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 500, 'TolFun', 1e-10, ...
                       'TolX', 1e-10, 'Diagnostics', 'off', 'Display', 'off');

    dnew = 0.1;
    B1 = (1 - 2.^(2*dnew-1))./(1-2*dnew);
    signew2 = va/(2*B1*(2*pi)^(1-2*dnew)*sum(2.^((2*dnew-1)*(Js))));

    it = 1;
    sigold2 = signew2/(1-abstol)+10;
    while (abs(1-signew2/sigold2) > abstol) && (it <= maxit)
        sigold2 = signew2;
        f = @(x) likelihood(x, S, Nj, [1, J], signew2);
        dnew = fminbnd(f, lb, ub, options);
        B1 = (1 - 2.^(2*dnew-1))./(1-2*dnew);
        signew2 = va/(2*B1*(2*pi)^(1-2*dnew)*sum(2.^((2*dnew-1)*(Js))));
        it = it + 1;
    end
    if it > 1
        d = dnew;
    else
        d = NaN;
    end
end

function d = profile_ml(S, Nj, J)
% Minimise the profile negative log-likelihood over d in [-0.5, 0.5].
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 500, 'TolFun', 1e-10, ...
                       'TolX', 1e-10, 'Diagnostics', 'off', 'Display', 'off');
    f = @(x) profile_nll(x, S, Nj, J);
    d = fminbnd(f, -0.5, 0.5, options);
end

function L = profile_nll(d, S, Nj, J)
    sdf = bfn_finsdf(d, [1, J], NaN);
    L   = sum(Nj) * log(sum(S ./ sdf)) + sum(Nj .* log(sdf));
end

function L = likelihood(d, S, Nj, range, sig)
    sdf = bfn_finsdf(d, range, NaN);
    Js  = (range(1):range(2))';
    L   = sum(Nj(Js).*log(sdf)) + sum(S(Js)./sdf)/sig;
end
