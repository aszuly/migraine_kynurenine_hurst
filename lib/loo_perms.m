function P = loo_perms(perm_idx, j)
% LOO_PERMS  Permutations of the n-1 remaining participants when participant
% j is left out, induced by the permutations of all n participants: j is
% deleted from every permutation and the larger indices are shifted down.
    [nperm, n] = size(perm_idx);
    T = perm_idx';
    T = reshape(T(T ~= j), n - 1, nperm)';
    P = T - (T > j);
end
