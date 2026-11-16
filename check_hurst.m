% check_hurst.m
%
% Recomputes the per-participant Hurst exponents from the primary parcel
% time series (data/denoised_variants.mat) and compares them with the
% deposited values in results/hurst_parcels.csv and
% results/hurst_networks.csv.
%
% The estimates are made with the same calls as in step3_main_analysis.m:
% bfn_mfin_ml of the nonfractal toolbox with its default options, and are
% also checked against the independent implementation
% lib/hurst_wavelet_ml.m. Network values are the means over the parcels of
% each network (data/parcels.csv), the whole-cortex value the mean over
% all 100 parcels.
%
% The script stops with an error if any recomputed value differs from
% the deposited one by more than TOL.
%
% Run from the repository folder:  matlab -batch check_hurst

clearvars
addpath('lib'); setup_paths();
TOL = 1e-12;

%% Inputs ------------------------------------------------------------------
V  = load(fullfile('data', 'denoised_variants.mat'));
PC = readtable(fullfile('data', 'parcels.csv'));
X  = V.(V.primary);
ids = cellstr(V.ids(:));
n  = size(X, 3);
assert(size(X, 2) == 100 && numel(ids) == n);
netw = PC.network_index(:)';
net_labels = {'VN','SMN','DAN','SVA','LIM','CON','DMN'};
nNet = numel(net_labels);
fprintf('check_hurst: %d participants, %d volumes, %d parcels; primary pipeline: %s\n', ...
        n, size(X, 1), size(X, 2), V.primary);

%% Hurst exponents (as in step 3) ------------------------------------------
H = zeros(n, 100); Hc = H;
for s = 1:n
    H(s, :)  = bfn_mfin_ml(X(:, :, s), 'verbose', 0, 'd_only', 1);
    Hc(s, :) = hurst_wavelet_ml(X(:, :, s), 'dwt', 'la8', 'variance');
end
HN = zeros(n, nNet);
for k = 1:nNet
    HN(:, k) = mean(H(:, netw == k), 2);
end
H_wb = mean(H, 2);

%% Comparison with the deposited values -------------------------------------
DP = readtable(fullfile('results', 'hurst_parcels.csv'));
DN = readtable(fullfile('results', 'hurst_networks.csv'));
assert(isequal(DP.id, ids) && isequal(DN.id, ids), ...
       'check_hurst: participant codes differ from those in results/.');
dP  = max(abs(H(:) - reshape(DP{:, compose('P%03d', 1:100)}, [], 1)));
dN  = max(abs(HN(:) - reshape(DN{:, net_labels}, [], 1)));
dWB = max(abs(H_wb - DN.whole_brain));
dC  = max(abs(H(:) - Hc(:)));

fprintf('toolbox estimate vs own implementation (lib/hurst_wavelet_ml.m): max |diff| = %g\n', dC);
fprintf('recomputed vs deposited, max |diff|:\n');
fprintf('  parcels (%d x 100, results/hurst_parcels.csv)      %g\n', n, dP);
fprintf('  networks (%d x %d, results/hurst_networks.csv)      %g\n', n, nNet, dN);
fprintf('  whole cortex (%d, results/hurst_networks.csv)       %g\n', n, dWB);
fprintf('network mean Hurst, mean (SD):');
for k = 1:nNet, fprintf('  %s %.3f (%.3f)', net_labels{k}, mean(HN(:, k)), std(HN(:, k))); end
fprintf('  whole cortex %.3f (%.3f)\n', mean(H_wb), std(H_wb));

if max([dP, dN, dWB, dC]) > TOL
    error('check_hurst:mismatch', ...
          'check_hurst: recomputed Hurst exponents differ from the deposited values (tolerance %g).', TOL);
end
fprintf('check_hurst: all Hurst exponents agree with the deposited values (tolerance %g).\n', TOL);
