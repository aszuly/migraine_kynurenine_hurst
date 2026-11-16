function [rho, p, rp] = ps_perm(x, y, Xd, perm_idx)
% PS_PERM  Partial Spearman correlation with a permutation p-value.
%
% The ranks of x are permuted across participants while each participant's
% y value and covariates stay together, which breaks the x-y association and
% keeps the association of y with the covariates. The two-sided p-value is
% (1 + number of permutations with |rho*| >= |rho|) / (nperm + 1)
% (Phipson & Smyth 2010).
    rx  = tiedrank(x);
    ry  = tiedrank(y);
    rho = resid_corr(rx, ry, Xd);
    rp  = resid_corr(rx(perm_idx'), ry, Xd);
    p   = (sum(abs(rp) >= abs(rho) - 1e-12) + 1) / (size(perm_idx, 1) + 1);
end
