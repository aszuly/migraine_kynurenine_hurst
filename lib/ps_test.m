function [rho, p, ci] = ps_test(x, y, Xd, R)
% PS_TEST  Partial Spearman correlation, permutation p-value and BCa interval.
%
%   [rho, p, ci] = ps_test(x, y, Xd, R)
%
% Xd holds the intercept and covariates; R comes from make_resamples (fields
% perm_idx, boot_idx, alpha). See ps_perm, ps_boot and bca.
    [rho, p] = ps_perm(x, y, Xd, R.perm_idx);
    if nargout > 2
        n  = numel(x);
        bs = ps_boot(x, y, Xd, R.boot_idx);
        jk = zeros(n, 1);
        for j = 1:n
            keep  = [1:j-1, j+1:n];
            jk(j) = ps_corr(x(keep), y(keep), Xd(keep, :));
        end
        ci = bca(rho, bs, jk, R.alpha);
    end
end
