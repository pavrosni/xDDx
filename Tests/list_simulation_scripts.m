function scriptTable = list_simulation_scripts()
%LIST_SIMULATION_SCRIPTS List scripts selectable by the regression runner.
%
%   list_simulation_scripts
%   scripts = list_simulation_scripts

supportFolder = fullfile(fileparts(mfilename('fullpath')), 'example_regression');
oldPath = path;
cleanup = onCleanup(@() path(oldPath));
addpath(supportFolder);

scriptTable = list_example_regression_scripts();
if nargout == 0
    disp(scriptTable);
    clear scriptTable;
end
end
