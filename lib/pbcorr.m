function r = pbcorr(x, y, beta)
% PBCORR  Percentage-bend correlation (Wilcox 1994), bending constant beta.
    if nargin < 3, beta = 0.2; end
    a = pbend(x(:), beta);
    b = pbend(y(:), beta);
    r = sum(a .* b) / sqrt(sum(a .^ 2) * sum(b .^ 2));
end

function a = pbend(x, beta)
    n     = numel(x);
    sx    = sort(abs(x - median(x)));
    omhat = sx(floor((1 - beta) * n));
    psi   = (x - median(x)) / omhat;
    i1    = sum(psi < -1);
    i2    = sum(psi > 1);
    xs    = sort(x);
    pbos  = (omhat * (i2 - i1) + sum(xs((i1 + 1):(n - i2)))) / (n - i1 - i2);
    a     = max(min((x - pbos) / omhat, 1), -1);
end
