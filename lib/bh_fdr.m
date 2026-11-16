function q = bh_fdr(p)
% BH_FDR  Benjamini-Hochberg adjusted p-values (q-values).
    p = p(:);
    m = numel(p);
    [ps, ord] = sort(p);
    qs = ps .* m ./ (1:m)';
    qs = flipud(min(cummin(flipud(qs)), 1));
    q = nan(m, 1);
    q(ord) = qs;
end
