function bs = ps_boot(x, y, Xd, boot_idx)
% PS_BOOT  Bootstrap replicates of the partial Spearman correlation.
% Ranks are recomputed within every bootstrap sample.
    B  = size(boot_idx, 1);
    RX = tiedrank(x(boot_idx'));
    RY = tiedrank(y(boot_idx'));
    bs = zeros(B, 1);
    for b = 1:B
        bs(b) = resid_corr(RX(:, b), RY(:, b), Xd(boot_idx(b, :), :));
    end
end
