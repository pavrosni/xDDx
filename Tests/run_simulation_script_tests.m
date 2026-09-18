function varargout = run_simulation_script_tests(varargin)
%RUN_SIMULATION_SCRIPT_TESTS Run selected reorganized simulation scripts.
%
% Examples:
%   run_simulation_script_tests('Scripts', 'fp_sf')
%   run_simulation_script_tests('Scripts', {'fp_sf', 'fp_transient'}, ...
%       'SimulationDevices', 'cpu')
%   run_simulation_script_tests('Scripts', 'all', ...
%       'SimulationDevices', 'both')
%   run_simulation_script_tests('Scripts', 'fp_sf', ...
%       'ReferenceProject', 'oldest')
%
% See LIST_SIMULATION_SCRIPTS for every valid script selection and
% RUN_EXAMPLE_REGRESSION_COMPARISON for all supported name-value options.

supportFolder = fullfile(fileparts(mfilename('fullpath')), 'example_regression');
oldPath = path;
cleanup = onCleanup(@() path(oldPath));
addpath(supportFolder);

if nargout == 0
    run_example_regression_comparison(varargin{:});
else
    [varargout{1:nargout}] = run_example_regression_comparison(varargin{:});
end
end
