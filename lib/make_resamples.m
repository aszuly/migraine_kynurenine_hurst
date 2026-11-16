function R = make_resamples(strata, seed, nperm, nboot)
% MAKE_RESAMPLES  Permutation and bootstrap indices shared by all tests.
%
%   R = make_resamples(strata, seed, nperm, nboot)
%
% strata is the logical aura indicator (n x 1). The random number generator
% is seeded here and the indices are drawn in a fixed order, so every script
% that calls this function with the same arguments uses identical resamples.
%
%   R.perm_idx      nperm x n   random permutations of 1..n
%   R.boot_idx      nboot x n   bootstrap samples drawn within aura strata
%   R.sub.perm_idx  nperm x m   permutations of the m participants without aura
%   R.sub.boot_idx  nboot x m   bootstrap samples of those m participants
%   R.alpha         0.05
    rng(seed, 'twister');
    strata = logical(strata(:));
    n = numel(strata);
    R.perm_idx = zeros(nperm, n);
    for k = 1:nperm
        R.perm_idx(k, :) = randperm(n);
    end
    R.boot_idx = zeros(nboot, n);
    for L = [false true]
        rows = find(strata == L);
        R.boot_idx(:, rows) = rows(randi(numel(rows), nboot, numel(rows)));
    end
    m = sum(~strata);
    R.sub.perm_idx = zeros(nperm, m);
    for k = 1:nperm
        R.sub.perm_idx(k, :) = randperm(m);
    end
    R.sub.boot_idx = randi(m, nboot, m);
    R.sub.alpha = 0.05;
    R.alpha = 0.05;
    R.seed = seed;
end
