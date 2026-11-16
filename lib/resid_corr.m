function r = resid_corr(a, b, Xd)
% RESID_CORR  Correlation of a and b after projecting out the columns of Xd.
% a may be a matrix (one variable per column); r is then a row vector.
    U  = orth(Xd);
    ea = a - U * (U' * a);
    eb = b - U * (U' * b);
    r  = (eb' * ea) ./ sqrt(sum(ea .^ 2, 1) * (eb' * eb));
end
