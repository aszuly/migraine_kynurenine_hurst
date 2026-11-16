% step0_export_inputs.m
%
% Exports the inputs of the analysis from the MATLAB workspace of the
% original analysis (../mi_hurst_kyn.mat) into plain files in which the
% participants are identified by codes. Run from this folder.
%
%   private/subject_key.csv           code <-> image folder.
%   data/participants.csv             code, aura status, age and the nine
%                                     plasma metabolite concentrations
%   data/parcels.csv                  Schaefer-100 parcel labels and Yeo
%                                     7-network assignment, from the
%                                     official label file
%   data/original_submission_ts.mat   parcel time series analysed in the
%                                     original submission (198 x 100 x 21),
%                                     kept only for comparison
%
% Participants are coded S01-S21 in the order of the original analysis.
% Aura status is taken from the SUBGROUP column of the clinical table
% (1 = migraine with aura, 2 = migraine without aura).

clearvars
SRC = fullfile('..', 'mi_hurst_kyn.mat');
LAB = fullfile('..', '..', 'toolboxes', 'schaefer', ...
               'Schaefer2018_100Parcels_7Networks_order.txt');
S = load(SRC, 'kyn', 'mwoa', 'mwoa_ts', 'netwidx');
kyn = S.kyn;
n = height(kyn);
assert(n == 21 && numel(S.mwoa) == n && size(S.mwoa_ts, 3) == n);
id = compose('S%02d', (1:n)');

%% Subject key (identifying) ---------------------------------------------
folder = cell(n, 1);
for i = 1:n
    folder{i} = ['mi_p_' erase(erase(S.mwoa(i).name, 'mwoa_'), '.txt')];
end
writetable(table(id, folder), fullfile('private', 'subject_key.csv'));

%% Participants -----------------------------------------------------------
vn  = kyn.Properties.VariableNames;
col = @(name) kyn{:, find(strcmpi(vn, name), 1)};
sg  = col('SUBGROUP(AURA)');
assert(all(ismember(sg, [1 2])));
aura = double(sg == 1);
P = table(id, aura, col('AGE'), col('TRP'), col('KYN'), col('KYNA'), col('ANA'), ...
          col('3HK'), col('XA'), col('3HANA'), col('PICA'), col('QUIN'), ...
    'VariableNames', {'id','aura','age','TRP','L_KYN','KYNA','ANA', ...
                      'HK3','XA','HANA3','PICA','QUINA'});
assert(~any(ismissing(P), 'all'));
writetable(P, fullfile('data', 'participants.csv'));
fprintf('participants: n = %d, %d with aura; age %d-%d years\n', ...
        n, sum(aura), min(P.age), max(P.age));

%% Parcel labels and networks ----------------------------------------------
L = readtable(LAB, 'FileType', 'text', 'ReadVariableNames', false, 'Delimiter', '\t');
label = string(L{:, 2});
assert(numel(label) == 100 && isequal(L{:, 1}, (1:100)'));
tok = regexp(label, '^7Networks_([LR]H)_([A-Za-z]+)_', 'tokens', 'once');
hemi = cellfun(@(c) c{1}, tok, 'UniformOutput', false);
yeo  = cellfun(@(c) c{2}, tok, 'UniformOutput', false);
yeo_names = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
abbrev    = {'VN','SMN','DAN','SVA','LIM','CON','DMN'};
[ok, network_index] = ismember(yeo, yeo_names);
assert(all(ok));
network = abbrev(network_index)';
network_original_submission = abbrev(S.netwidx(:))';
parcel = (1:100)';
T = table(parcel, cellstr(label), hemi, network, network_index, network_original_submission, ...
    'VariableNames', {'parcel','label','hemisphere','network','network_index', ...
                      'network_original_submission'});
writetable(T, fullfile('data', 'parcels.csv'));
changed = find(~strcmp(network, network_original_submission));
fprintf('parcels per network (official labels): %s\n', mat2str(accumarray(network_index, 1)'));
fprintf('parcels whose network differs from the original submission: %d (%s)\n', ...
        numel(changed), mat2str(changed'));

%% Series analysed in the original submission -------------------------------
ts = S.mwoa_ts;
save(fullfile('data', 'original_submission_ts.mat'), 'ts');
fprintf('done\n');
