% step3_main_analysis.m
%
% Main analysis, on the primary parcel time series
% (data/denoised_variants.mat). The permutation and bootstrap resamples are
% drawn once (lib/make_resamples.m), and step4_sensitivity.m uses the same
% ones.
%
% Hurst exponent: nonfractal toolbox, bfn_mfin_ml with default options
%   (DWT, LA8, first 128 of 198 volumes, levels 1-4, maximum likelihood).
% Association:    partial Spearman correlation, covariates age and aura.
% p-values:       20,000 permutations of the metabolite across participants.
% Intervals:      BCa bootstrap, 10,000 resamples stratified by aura status.
% Multiplicity:   Benjamini-Hochberg within the primary family (L-KYN x 7
%                 networks), the secondary family (9 other measures x 7
%                 networks) and across all 70 tests.
%
% Outputs in results/ (and formatted supplementary tables in tables/):
%   hurst_parcels.csv, hurst_networks.csv, tests_70.csv, wholebrain.csv,
%   network_diff.csv, smn_lkyn_robustness.csv, mwoa_only.csv, loo.csv,
%   percentage_bend.csv, motion_tests.csv, ratios.csv,
%   network_label_check.csv, fig_subjects.csv, log_step3.txt

clearvars
addpath('lib'); setup_paths();
SEED = 20261105; NPERM = 20000; NBOOT = 10000;
if exist(fullfile('results', 'log_step3.txt'), 'file'), delete(fullfile('results', 'log_step3.txt')); end
diary(fullfile('results', 'log_step3.txt'));

%% Inputs ------------------------------------------------------------------
P  = readtable(fullfile('data', 'participants.csv'));
PC = readtable(fullfile('data', 'parcels.csv'));
MO = readtable(fullfile('data', 'motion.csv'));
V  = load(fullfile('data', 'denoised_variants.mat'));
X  = V.(V.primary);
n  = height(P);
assert(isequal(P.id, MO.id) && size(X, 3) == n && size(X, 2) == 100);
aura = logical(P.aura);
Xd   = [ones(n, 1), P.age, double(aura)];
netw = PC.network_index(:)';
net_labels = {'VN','SMN','DAN','SVA','LIM','CON','DMN'};
nNet = numel(net_labels);
iS   = find(strcmp(net_labels, 'SMN'));

met_cols  = {'TRP','L_KYN','KYNA','ANA','HK3','XA','HANA3','PICA','QUINA'};
met_names = {'TRP','L-KYN','KYNA','ANA','3-HK','XA','3-HANA','PICA','QUINA','QUINA/KYNA'};
M = [P{:, met_cols}, P.QUINA ./ P.KYNA];
nMet = numel(met_names);
LK = P.L_KYN;
fprintf('n = %d (%d with aura, %d without); primary pipeline: %s\n', n, sum(aura), sum(~aura), V.primary);

R = make_resamples(aura, SEED, NPERM, NBOOT);

%% Hurst exponents ---------------------------------------------------------
H = zeros(n, 100); Hc = H;
for s = 1:n
    H(s, :)  = bfn_mfin_ml(X(:, :, s), 'verbose', 0, 'd_only', 1);
    Hc(s, :) = hurst_wavelet_ml(X(:, :, s), 'dwt', 'la8', 'variance');
end
fprintf('toolbox estimate vs own implementation: max |diff| = %g\n', max(abs(H(:) - Hc(:))));
HN = zeros(n, nNet);
for k = 1:nNet
    HN(:, k) = mean(H(:, netw == k), 2);
end
H_wb = mean(H, 2);
fprintf('parcel estimates at a bound (H < 0.001 or > 0.999): %.1f%%\n', 100 * mean(H(:) < 1e-3 | H(:) > 1 - 1e-3));
fprintf('network mean Hurst, mean (SD):');
for k = 1:nNet, fprintf('  %s %.3f (%.3f)', net_labels{k}, mean(HN(:, k)), std(HN(:, k))); end
fprintf('  whole brain %.3f (%.3f)\n', mean(H_wb), std(H_wb));

T = array2table(H, 'VariableNames', compose('P%03d', 1:100));
writetable([P(:, 'id'), T], fullfile('results', 'hurst_parcels.csv'));
T = array2table([HN, H_wb], 'VariableNames', [net_labels, {'whole_brain'}]);
writetable([P(:, 'id'), T], fullfile('results', 'hurst_networks.csv'));
writetable(table(P.id, double(aura), LK, HN(:, iS), 'VariableNames', {'id','aura','L_KYN','SMN_Hurst'}), ...
           fullfile('results', 'fig_subjects.csv'));

%% 70 network x measure tests ----------------------------------------------
rows = cell(nNet * nMet, 6); c = 0;
for m = 1:nMet
    for k = 1:nNet
        c = c + 1;
        [rho, p, ci] = ps_test(M(:, m), HN(:, k), Xd, R);
        rows(c, :) = {net_labels{k}, met_names{m}, rho, ci(1), ci(2), p};
    end
end
T70 = cell2table(rows, 'VariableNames', {'Network','Measure','rho','CI_lo','CI_hi','p_perm'});
prim = strcmp(T70.Measure, 'L-KYN');
T70.family = repmat({'secondary'}, height(T70), 1);
T70.family(prim) = {'primary'};
T70.q_family = nan(height(T70), 1);
T70.q_family(prim)  = bh_fdr(T70.p_perm(prim));
T70.q_family(~prim) = bh_fdr(T70.p_perm(~prim));
T70.q_all70 = bh_fdr(T70.p_perm);
writetable(T70, fullfile('results', 'tests_70.csv'));
writetable(T70, fullfile('tables', 'Table_S1_network_by_measure.csv'));
fprintf('\n=== L-KYN x network (primary family) ===\n');
disp(T70(prim, {'Network','rho','CI_lo','CI_hi','p_perm','q_family','q_all70'}));
fprintf('Secondary family: smallest p = %.4f, smallest q = %.4f; tests with p < 0.05: %d of 63\n', ...
        min(T70.p_perm(~prim)), min(T70.q_family(~prim)), sum(T70.p_perm(~prim) < 0.05));
[~, ord] = sort(T70.p_perm);
fprintf('Five smallest p-values of all 70:\n');
disp(T70(ord(1:5), {'Network','Measure','rho','p_perm','q_family','q_all70'}));

%% Whole-brain mean Hurst ----------------------------------------------------
rows = cell(nMet, 5);
for m = 1:nMet
    [rho, p, ci] = ps_test(M(:, m), H_wb, Xd, R);
    rows(m, :) = {met_names{m}, rho, ci(1), ci(2), p};
end
WB = cell2table(rows, 'VariableNames', {'Measure','rho','CI_lo','CI_hi','p_perm'});
writetable(WB, fullfile('results', 'wholebrain.csv'));
fprintf('\n=== Whole-brain mean Hurst ===\n'); disp(WB);

%% SMN against the other networks, same participants ------------------------
Yall = [HN, mean(H(:, netw ~= iS), 2)];
ylab = [net_labels, {'all other parcels'}];
BS = zeros(NBOOT, numel(ylab)); JK = zeros(n, numel(ylab)); OB = zeros(1, numel(ylab));
for k = 1:numel(ylab)
    OB(k) = ps_corr(LK, Yall(:, k), Xd);
    BS(:, k) = ps_boot(LK, Yall(:, k), Xd, R.boot_idx);
    for j = 1:n
        keep = [1:j-1, j+1:n];
        JK(j, k) = ps_corr(LK(keep), Yall(keep, k), Xd(keep, :));
    end
end
rows = {};
for k = setdiff(1:numel(ylab), iS)
    d_bs = BS(:, iS) - BS(:, k);
    ci = bca(OB(iS) - OB(k), d_bs, JK(:, iS) - JK(:, k), R.alpha);
    p_boot = min(1, 2 * min(mean(d_bs <= 0), mean(d_bs >= 0)));
    rows(end + 1, :) = {sprintf('SMN - %s', ylab{k}), OB(k), OB(iS) - OB(k), ci(1), ci(2), p_boot}; %#ok<SAGROW>
end
ND = cell2table(rows, 'VariableNames', {'Contrast','rho_other','diff_rho','CI_lo','CI_hi','p_boot'});
writetable(ND, fullfile('results', 'network_diff.csv'));
fprintf('\n=== Difference in the L-KYN partial rho: SMN minus other network ===\n'); disp(ND);

%% Migraine without aura only -------------------------------------------------
w = ~aura; Xw = [ones(sum(w), 1), P.age(w)];
rows = {};
for k = 1:nNet + 1
    if k <= nNet, y = HN(w, k); lab = net_labels{k}; else, y = H_wb(w); lab = 'whole brain'; end
    [rho, p, ci] = ps_test(LK(w), y, Xw, R.sub);
    rows(end + 1, :) = {lab, sum(w), rho, ci(1), ci(2), p}; %#ok<SAGROW>
end
MW = cell2table(rows, 'VariableNames', {'Network','n','rho','CI_lo','CI_hi','p_perm'});
writetable(MW, fullfile('results', 'mwoa_only.csv'));
fprintf('\n=== L-KYN x network, migraine without aura only (covariate: age) ===\n'); disp(MW);

%% Leave-one-out ---------------------------------------------------------------
rows = cell(n * (nNet + 1), 5); c = 0;
for k = 1:nNet + 1
    if k <= nNet, y = HN(:, k); lab = net_labels{k}; else, y = H_wb; lab = 'whole brain'; end
    for j = 1:n
        keep = [1:j-1, j+1:n];
        [rho, p] = ps_perm(LK(keep), y(keep), Xd(keep, :), loo_perms(R.perm_idx, j));
        c = c + 1;
        rows(c, :) = {lab, 'L-KYN', P.id{j}, rho, p};
    end
end
LOO = cell2table(rows, 'VariableNames', {'Network','Measure','left_out','rho','p_perm'});
writetable(LOO, fullfile('results', 'loo.csv'));
writetable(LOO, fullfile('tables', 'Table_S2_leave_one_out.csv'));
fprintf('\n=== Leave-one-out, L-KYN x network: rho range (median), largest p ===\n');
for k = 1:nNet + 1
    if k <= nNet, lab = net_labels{k}; else, lab = 'whole brain'; end
    q = LOO(strcmp(LOO.Network, lab), :);
    fprintf('  %-12s %.2f to %.2f (median %.2f), p %.4f to %.4f\n', lab, min(q.rho), max(q.rho), ...
            median(q.rho), min(q.p_perm), max(q.p_perm));
end

%% Percentage-bend correlation on covariate-adjusted values --------------------
res = @(v) v - Xd * (Xd \ v);
rows = {};
for k = 1:nNet + 1
    if k <= nNet, y = HN(:, k); lab = net_labels{k}; else, y = H_wb; lab = 'whole brain'; end
    ey = res(y);
    r  = pbcorr(res(LK), ey, 0.2);
    rp = zeros(NPERM, 1);
    for b = 1:NPERM
        rp(b) = pbcorr(res(LK(R.perm_idx(b, :))), ey, 0.2);
    end
    p_perm = (sum(abs(rp) >= abs(r) - 1e-12) + 1) / (NPERM + 1);
    tstat  = r * sqrt((n - 4) / (1 - r ^ 2));
    rows(end + 1, :) = {lab, r, p_perm, 2 * tcdf(-abs(tstat), n - 4)}; %#ok<SAGROW>
end
PB = cell2table(rows, 'VariableNames', {'Network','r_pb','p_perm','p_t_df_n_minus_4'});
writetable(PB, fullfile('results', 'percentage_bend.csv'));
fprintf('\n=== Percentage-bend correlation (beta = 0.2), L-KYN x network, adjusted for age and aura ===\n'); disp(PB);

%% Head motion -------------------------------------------------------------------
FDJ = MO.FD_jenkinson_mean; FDP = MO.FD_power_mean; AF = MO.aroma_motion_ICs ./ MO.ICs_total;
rows = {};
mot = {'mean FD (Jenkinson)', FDJ; 'mean FD (Power)', FDP; 'AROMA motion-component share', AF};
for i = 1:size(mot, 1)
    tgt = {'L-KYN', LK; 'SMN Hurst', HN(:, iS); 'whole-brain Hurst', H_wb};
    for t = 1:size(tgt, 1)
        [rho, p, ci] = ps_test(mot{i, 2}, tgt{t, 2}, Xd, R);
        rows(end + 1, :) = {mot{i, 1}, tgt{t, 1}, 'age, aura', rho, ci(1), ci(2), p}; %#ok<SAGROW>
    end
end
cov = {'age, aura, mean FD (Jenkinson)', [Xd, FDJ]; 'age, aura, mean FD (Power)', [Xd, FDP]};
for i = 1:size(cov, 1)
    for k = 1:nNet + 1
        if k <= nNet, y = HN(:, k); lab = [net_labels{k} ' Hurst']; else, y = H_wb; lab = 'whole-brain Hurst'; end
        [rho, p, ci] = ps_test(LK, y, cov{i, 2}, R);
        rows(end + 1, :) = {'L-KYN', lab, cov{i, 1}, rho, ci(1), ci(2), p}; %#ok<SAGROW>
    end
end
MT = cell2table(rows, 'VariableNames', {'x','y','covariates','rho','CI_lo','CI_hi','p_perm'});
writetable(MT, fullfile('results', 'motion_tests.csv'));
writetable(MT, fullfile('tables', 'Table_S6_head_motion.csv'));
fprintf('\n=== Head motion ===\n'); disp(MT);

%% Other covariate sets for SMN x L-KYN ---------------------------------------------
rows = {};
sets = {'age, aura (main model)', Xd; 'none', ones(n, 1); 'age', Xd(:, 1:2); 'aura', Xd(:, [1 3])};
for i = 1:size(sets, 1)
    [rho, p, ci] = ps_test(LK, HN(:, iS), sets{i, 2}, R);
    rows(end + 1, :) = {sets{i, 1}, n, rho, ci(1), ci(2), p}; %#ok<SAGROW>
end
q = MT(strcmp(MT.x, 'L-KYN') & strcmp(MT.y, 'SMN Hurst'), :);
for i = 1:height(q)
    rows(end + 1, :) = {q.covariates{i}, n, q.rho(i), q.CI_lo(i), q.CI_hi(i), q.p_perm(i)}; %#ok<SAGROW>
end
rows(end + 1, :) = {'age; migraine without aura only', MW.n(iS), MW.rho(iS), MW.CI_lo(iS), MW.CI_hi(iS), MW.p_perm(iS)};
RB = cell2table(rows, 'VariableNames', {'covariates','n','rho','CI_lo','CI_hi','p_perm'});
writetable(RB, fullfile('results', 'smn_lkyn_robustness.csv'));
fprintf('\n=== SMN x L-KYN under other covariate sets ===\n'); disp(RB);

%% Metabolite ratios -------------------------------------------------------------------
ratios = {'KYN/TRP', P.L_KYN ./ P.TRP; 'KYNA/KYN', P.KYNA ./ P.L_KYN; ...
          '3-HK/KYN', P.HK3 ./ P.L_KYN; 'QUINA/KYNA', P.QUINA ./ P.KYNA};
rows = {};
for i = 1:size(ratios, 1)
    for k = 1:nNet
        [rho, p, ci] = ps_test(ratios{i, 2}, HN(:, k), Xd, R);
        rows(end + 1, :) = {ratios{i, 1}, net_labels{k}, rho, ci(1), ci(2), p}; %#ok<SAGROW>
    end
end
RT = cell2table(rows, 'VariableNames', {'Ratio','Network','rho','CI_lo','CI_hi','p_perm'});
RT.q_28 = bh_fdr(RT.p_perm);
writetable(RT, fullfile('results', 'ratios.csv'));
writetable(RT, fullfile('tables', 'Table_S7_metabolite_ratios.csv'));
fprintf('\n=== Metabolite ratios x network (28 tests) ===\n'); disp(RT);

%% Effect of correcting the network labels -------------------------------------------------
[~, old_idx] = ismember(PC.network_original_submission, net_labels);
rows = {};
for k = 1:nNet
    [r_new, p_new] = ps_perm(LK, HN(:, k), Xd, R.perm_idx);
    [r_old, p_old] = ps_perm(LK, mean(H(:, old_idx == k), 2), Xd, R.perm_idx);
    rows(end + 1, :) = {net_labels{k}, sum(netw == k), r_new, p_new, sum(old_idx == k), r_old, p_old}; %#ok<SAGROW>
end
NL = cell2table(rows, 'VariableNames', {'Network','parcels_official','rho_official','p_official', ...
                                        'parcels_original','rho_original_labels','p_original_labels'});
writetable(NL, fullfile('results', 'network_label_check.csv'));
fprintf('\n=== L-KYN x network with official vs originally used parcel labels ===\n'); disp(NL);

save(fullfile('results', 'step3_workspace.mat'), 'H', 'HN', 'H_wb', 'SEED', 'NPERM', 'NBOOT');
diary off
