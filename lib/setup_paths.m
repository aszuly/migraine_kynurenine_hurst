function setup_paths()
% SETUP_PATHS  Add the Hurst toolboxes (nonfractal, WMTSA and the compiled
% WMTSA MEX files) to the MATLAB path; see toolboxes/README.md.
%
% The toolbox folder is taken from the environment variable HURST_TOOLBOXES
% if it is set, otherwise from toolboxes/ in this repository.
    here = fileparts(mfilename('fullpath'));
    TB = getenv('HURST_TOOLBOXES');
    if isempty(TB)
        TB = fullfile(here, '..', 'toolboxes');
    end
    need = {'nonfractal', 'wmtsa', 'wmtsa_mex'};
    for k = 1:numel(need)
        if ~isfolder(fullfile(TB, need{k}))
            error('setup_paths:missing', ...
                  'Toolbox folder %s not found in %s; see toolboxes/README.md.', ...
                  need{k}, TB);
        end
    end
    addpath(genpath(fullfile(TB, 'nonfractal')));
    addpath(genpath(fullfile(TB, 'wmtsa')));
    addpath(fullfile(TB, 'wmtsa_mex'), '-begin');
end
