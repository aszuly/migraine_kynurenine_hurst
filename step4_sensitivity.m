% step4_sensitivity.m
%
% Sensitivity analyses for the L-KYN x Hurst association. Uses the same
% seed, permutations and bootstrap samples as step3_main_analysis.m.
%
%   A  Wavelet estimator settings (primary time series)
%   B  Denoising variants (default estimator)
%   C  Position of the 128-volume window (default and shift-invariant estimator)
%   D  Additional checks: AROMA motion-component share as a covariate;
%      parcels partly outside the brain mask excluded; participants with
%      the most head motion excluded
%
% Outputs in results/: wavelet_settings.csv, wavelet_settings_summary.csv,
% denoising_variants.csv, denoising_variants_summary.csv, windows.csv,
% extra_checks.csv, log_step4.txt; supplementary tables S3-S5 and S8 in tables/.

clearvars
addpath('lib'); setup_paths();
SEED = 20261105; NPERM = 20000; NBOOT = 10000;
if exist(fullfile('results', 'log_step4.txt'), 'file'), delete(fullfile('results', 'log_step4.txt')); end
diary(fullfile('results', 'log_step4.txt'));

P  = readtable(fullfile('data', 'participants.csv'));
PC = readtable(fullfile('data', 'parcels.csv'));
MO = readtable(fullfile('data', 'motion.csv'));
CV = readtable(fullfile('data', 'parcel_coverage.csv'));
V  = load(fullfile('data', 'denoised_variants.mat'));
S3 = load(fullfile('results', 'step3_workspace.mat'));
X  = V.(V.primary);
n  = height(P);
assert(isequal(P.id, MO.id) && isequal(P.id, CV.id));
aura = logical(P.aura);
Xd   = [ones(n, 1), P.age, double(aura)];
netw = PC.network_index(:)';
net_labels = {'VN','SMN','DAN','SVA','LIM','CON','DMN'};
targets = [net_labels, {'whole brain'}];
LK = P.L_KYN;
R  = make_resamples(aura, SEED, NPERM, NBOOT);

%% A. Wavelet estimator settings ---------------------------------------------
tb = @(X, varargin) bfn_mfin_ml(X, 'verbose', 0, 'd_only', 1, varargin{:});
settings = {
    'DWT, LA8, levels 1-4, volumes 1-128 (main analysis)',    1:198,  @(X) tb(X)
    'DWT, LA8, levels 1-3, volumes 1-128',                    1:128,  @(X) tb(X, 'range', [1 3])
    'DWT, LA8, levels 2-4, volumes 1-128',                    1:128,  @(X) tb(X, 'range', [2 4])
    'DWT, Haar, levels 1-4, volumes 1-128',                   1:128,  @(X) tb(X, 'filter', 'haar', 'range', [1 4])
    'DWT, D4, levels 1-4, volumes 1-128',                     1:128,  @(X) tb(X, 'filter', 'd4', 'range', [1 4])
    'DWT, C6, levels 1-4, volumes 1-128',                     1:128,  @(X) tb(X, 'filter', 'c6', 'range', [1 4])
    'DWT, LA8, levels 1-4, volumes 71-198',                   71:198, @(X) tb(X)
    'DWT, LA8, levels 1-4, volumes 1-128, profile likelihood', 1:128, @(X) hurst_wavelet_ml(X, 'dwt_eff', 'la8', 'profile', 4)
    'MODWT, LA8, levels 1-4, volumes 1-128',                  1:128,  @(X) hurst_wavelet_ml(X, 'modwt', 'la8', 'profile', 4)
    'MODWT, LA8, levels 1-4, all 198 volumes',                1:198,  @(X) hurst_wavelet_ml(X, 'modwt', 'la8', 'profile', 4)
    'MODWT, Haar, levels 1-4, all 198 volumes',               1:198,  @(X) hurst_wavelet_ml(X, 'modwt', 'haar', 'profile', 4)
    'MODWT, D4, levels 1-4, all 198 volumes',                 1:198,  @(X) hurst_wavelet_ml(X, 'modwt', 'd4', 'profile', 4)
    };
long = {}; summ = {};
for i = 1:size(settings, 1)
    H = zeros(n, 100);
    for s = 1:n
        H(s, :) = settings{i, 3}(X(settings{i, 2}, :, s));
    end
    if i == 1
        assert(isequal(H, S3.H), 'main-analysis estimates differ from step 3');
    end
    [rows, srow] = eval_hurst(H, netw, LK, Xd, R, targets, S3.HN(:, 2));
    long = [long; [repmat(settings(i, 1), size(rows, 1), 1), rows]]; %#ok<AGROW>
    summ = [summ; [settings(i, 1), srow]]; %#ok<AGROW>
    fprintf('A %-58s SMN rho %.3f [%.2f, %.2f] p %.4f | whole brain %.3f p %.4f | at bounds %.1f%% / %.1f%%\n', ...
            settings{i, 1}, rows{2, 2}, rows{2, 3}, rows{2, 4}, rows{2, 5}, rows{8, 2}, rows{8, 5}, srow{2}, srow{3});
end
WS = finish_tables(long, summ, 'Setting');
writetable(WS.long, fullfile('results', 'wavelet_settings.csv'));
writetable(WS.summary, fullfile('results', 'wavelet_settings_summary.csv'));
writetable(WS.wide, fullfile('tables', 'Table_S3_wavelet_settings.csv'));

%% B. Denoising variants ---------------------------------------------------------
vlab = {
    'f24_hp',                 'Friston-24 regression, then high-pass (main analysis)'
    'hp_then_f24',            'High-pass (data and regressors), then Friston-24 regression'
    'f24deriv_hp',            '24 parameters as 6 + derivatives + squares, then high-pass'
    'f24_wmcsf_hp',           'Friston-24 + white-matter and CSF signals, then high-pass'
    'mp6_hp',                 '6 realignment parameters, then high-pass'
    'mp6_wmcsf_hp',           '6 realignment parameters + white-matter and CSF signals, then high-pass'
    'mcf_hp',                 'No nuisance regression, high-pass only'
    'aroma_nonaggr_wmcsf_hp', 'ICA-AROMA (non-aggressive) + white-matter and CSF signals, then high-pass'
    'aroma_aggr_hp',          'ICA-AROMA motion components regressed out in full, then high-pass'
    'original_submission',    'Time series analysed in the original submission'
    };
long = {}; summ = {};
for i = 1:size(vlab, 1)
    Xv = V.(vlab{i, 1});
    H = zeros(n, 100);
    for s = 1:n
        H(s, :) = tb(Xv(:, :, s));
    end
    if i == 1
        assert(isequal(H, S3.H), 'main-analysis estimates differ from step 3');
    end
    [rows, srow] = eval_hurst(H, netw, LK, Xd, R, targets, S3.HN(:, 2));
    long = [long; [repmat(vlab(i, 2), size(rows, 1), 1), rows]]; %#ok<AGROW>
    summ = [summ; [vlab(i, 2), srow]]; %#ok<AGROW>
    fprintf('B %-74s SMN rho %.3f [%.2f, %.2f] p %.4f | whole brain %.3f p %.4f | mean H %.2f\n', ...
            vlab{i, 2}, rows{2, 2}, rows{2, 3}, rows{2, 4}, rows{2, 5}, rows{8, 2}, rows{8, 5}, srow{1});
end
DV = finish_tables(long, summ, 'Denoising');
writetable(DV.long, fullfile('results', 'denoising_variants.csv'));
writetable(DV.summary, fullfile('results', 'denoising_variants_summary.csv'));
writetable(DV.wide, fullfile('tables', 'Table_S4_denoising_variants.csv'));

%% C. Position of the 128-volume window ---------------------------------------------
% Default estimator (DWT; depends on where the window starts) and the
% shift-invariant MODWT estimator with the same filter and levels.
west = {'DWT, LA8, levels 1-4 (main analysis)', @(X) tb(X)
        'MODWT, LA8, levels 1-4',               @(X) hurst_wavelet_ml(X, 'modwt', 'la8', 'profile', 4)};
rows = {};
for e = 1:size(west, 1)
    Y1 = [];
    for st = 1:10:71
        H = zeros(n, 100);
        for s = 1:n
            H(s, :) = west{e, 2}(X(st:st + 127, :, s));
        end
        Y = net_means(H, netw);
        if st == 1, Y1 = Y; end
        for k = 1:numel(targets)
            [rho, p] = ps_perm(LK, Y(:, k), Xd, R.perm_idx);
            rows(end + 1, :) = {west{e, 1}, st, st + 127, targets{k}, rho, p, ...
                                corr(Y(:, k), Y1(:, k), 'type', 'Spearman')}; %#ok<SAGROW>
        end
    end
end
WN = cell2table(rows, 'VariableNames', {'Estimator','first_volume','last_volume','Network','rho','p_perm', ...
                                        'rho_H_with_first_window'});
writetable(WN, fullfile('results', 'windows.csv'));
writetable(WN(ismember(WN.Network, {'SMN','whole brain'}), :), fullfile('tables', 'Table_S5_window_position.csv'));
fprintf('\nC. 128-volume windows, L-KYN partial rho (p); Spearman of network Hurst with the first window\n');
for e = 1:size(west, 1)
    fprintf(' %s\n', west{e, 1});
    for st = 1:10:71
        q = WN(WN.first_volume == st & strcmp(WN.Estimator, west{e, 1}), :);
        fprintf('  volumes %3d-%3d:', st, st + 127);
        for k = [2 8], fprintf('  %s %5.2f (%.4f), with first window %.2f;', q.Network{k}, q.rho(k), q.p_perm(k), q.rho_H_with_first_window(k)); end
        fprintf('\n');
    end
end

%% D. Additional checks -----------------------------------------------------------------
AF  = MO.aroma_motion_ICs ./ MO.ICs_total;
cov = CV{:, 2:end};
low = any(cov < 0.5, 1);
fprintf('\nD. parcels with < 50%% inside the brain mask in at least one participant: %s (networks %s)\n', ...
        mat2str(find(low)), strjoin(net_labels(netw(low)), ', '));
Yall = net_means(S3.H, netw);
Hx = S3.H; Hx(:, low) = NaN;
Ycov = zeros(n, 8);
for k = 1:7, Ycov(:, k) = mean(Hx(:, netw == k), 2, 'omitnan'); end
Ycov(:, 8) = mean(Hx, 2, 'omitnan');
rows = {};
for k = 1:numel(targets)
    [rho, p, ci] = ps_test(LK, Yall(:, k), [Xd, AF], R);
    rows(end + 1, :) = {'AROMA motion-component share added as covariate', targets{k}, rho, ci(1), ci(2), p}; %#ok<SAGROW>
end
for k = 1:numel(targets)
    [rho, p, ci] = ps_test(LK, Ycov(:, k), Xd, R);
    rows(end + 1, :) = {'parcels < 50% inside the brain mask excluded', targets{k}, rho, ci(1), ci(2), p}; %#ok<SAGROW>
end
% Participants with more than 5% of volumes above FD (Power) = 0.5 mm left
% out. The subset needs its own resamples. They are drawn last, from seed
% SEED + 1, so no other result depends on them.
hi = MO.pct_vols_FDpower_gt_0p5 > 5;
fprintf('participants with > 5%% of volumes above FD (Power) 0.5 mm: %s (mean FD, Jenkinson: %s mm)\n', ...
        strjoin(P.id(hi)', ', '), strjoin(compose('%.2f', MO.FD_jenkinson_mean(hi))', ', '));
Rk = make_resamples(aura(~hi), SEED + 1, NPERM, NBOOT);
lab = sprintf('participants with > 5%% of volumes above FD 0.5 mm excluded (n = %d)', sum(~hi));
for k = 1:numel(targets)
    [rho, p, ci] = ps_test(LK(~hi), Yall(~hi, k), Xd(~hi, :), Rk);
    rows(end + 1, :) = {lab, targets{k}, rho, ci(1), ci(2), p}; %#ok<SAGROW>
end
EX = cell2table(rows, 'VariableNames', {'Check','Network','rho','CI_lo','CI_hi','p_perm'});
writetable(EX, fullfile('results', 'extra_checks.csv'));
writetable(EX, fullfile('tables', 'Table_S8_additional_checks.csv'));
disp(EX);
diary off

%% Local functions -----------------------------------------------------------------------
function Y = net_means(H, netw)
    nNet = max(netw);
    Y = zeros(size(H, 1), nNet + 1);
    for k = 1:nNet
        Y(:, k) = mean(H(:, netw == k), 2);
    end
    Y(:, nNet + 1) = mean(H, 2);
end

function [rows, srow] = eval_hurst(H, netw, LK, Xd, R, targets, smn_main)
    Y = net_means(H, netw);
    rows = cell(numel(targets), 5);
    for k = 1:numel(targets)
        [rho, p, ci] = ps_test(LK, Y(:, k), Xd, R);
        rows(k, :) = {targets{k}, rho, ci(1), ci(2), p};
    end
    srow = {mean(H(:)), 100 * mean(H(:) < 1e-3), 100 * mean(H(:) > 1 - 1e-3), ...
            corr(Y(:, 2), smn_main, 'type', 'Spearman')};
end

function out = finish_tables(long, summ, name)
    L = cell2table(long, 'VariableNames', {name, 'Network', 'rho', 'CI_lo', 'CI_hi', 'p_perm'});
    L.q_7networks = nan(height(L), 1);
    labs = unique(L.(name), 'stable');
    for i = 1:numel(labs)
        sel = strcmp(L.(name), labs{i}) & ~strcmp(L.Network, 'whole brain');
        L.q_7networks(sel) = bh_fdr(L.p_perm(sel));
    end
    S = cell2table(summ, 'VariableNames', {name, 'mean_H', 'pct_at_lower_bound', 'pct_at_upper_bound', ...
                                           'rho_SMN_H_with_main_analysis'});
    W = S;
    nets = unique(L.Network, 'stable');
    for k = 1:numel(nets)
        q = L(strcmp(L.Network, nets{k}), :);
        tag = strrep(nets{k}, ' ', '_');
        W.([tag '_rho']) = q.rho; W.([tag '_CI_lo']) = q.CI_lo; W.([tag '_CI_hi']) = q.CI_hi;
        W.([tag '_p']) = q.p_perm;
        if ~strcmp(nets{k}, 'whole brain'), W.([tag '_q7']) = q.q_7networks; end
    end
    out = struct('long', L, 'summary', S, 'wide', W);
end
