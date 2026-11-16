function r = ps_corr(x, y, Xd)
% PS_CORR  Partial Spearman correlation of x and y given Xd.
% Xd contains the intercept and the covariates. x and y are ranked, the
% column space of Xd is projected out of both rank vectors, and the Pearson
% correlation of the residuals is returned.
    r = resid_corr(tiedrank(x), tiedrank(y), Xd);
end
