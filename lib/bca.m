function ci = bca(obs, bs, jk, alpha)
% BCA  Bias-corrected and accelerated bootstrap interval (Efron 1987).
% obs: observed statistic; bs: bootstrap replicates; jk: jackknife values.
    bs = bs(~isnan(bs));
    z0 = norminv(max(min(mean(bs < obs), 1 - eps), eps));
    jm = mean(jk);
    a  = sum((jm - jk) .^ 3) / (6 * (sum((jm - jk) .^ 2)) ^ 1.5 + eps);
    zl = norminv(alpha / 2);
    zu = norminv(1 - alpha / 2);
    q1 = normcdf(z0 + (z0 + zl) / (1 - a * (z0 + zl)));
    q2 = normcdf(z0 + (z0 + zu) / (1 - a * (z0 + zu)));
    ci = quantile(bs, [q1 q2]);
end
