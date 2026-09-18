function [results, reportPaths, hostInfo] = ...
    run_xddx_backend_comparison(options)
%RUN_XDDX_BACKEND_COMPARISON Run the default simulator case on all backends.
%   The numerical setup comes directly from xDDx_simulator.m. Only backend,
%   launch mode, GUI, and plotting settings are overridden.

if nargin < 1
    options = xddx_backend_options_from_environment();
end

testFolder = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(testFolder));
simulatorFolder = fullfile(projectRoot, 'simulation_toolbox', ...
    'heterogeneous_simulator');
simulatorScript = fullfile(simulatorFolder, 'xDDx_simulator.m');

addpath(fullfile(projectRoot, 'xDDx_lib'));
addpath(testFolder);
cases = xddx_backend_matrix(options);
hostInfo = xddx_backend_host_info(projectRoot);
results = repmat(empty_result(), 1, numel(cases));
runInfo = create_run_info();
reportPaths = struct('Csv', '', 'Mat', '');

oldXddxDocker = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
oldKwaveDocker = getenv('KWAVE_USE_DOCKER_ON_WINDOWS');
environmentCleanup = onCleanup(@() restore_environment( ...
    oldXddxDocker, oldKwaveDocker));

reference = [];
referenceCase = '';
for caseIndex = 1:numel(cases)
    currentCase = cases(caseIndex);
    results(caseIndex) = result_for_case(currentCase);
    [available, skipReason] = xddx_backend_case_availability( ...
        currentCase, hostInfo, projectRoot);
    if ~available
        results(caseIndex).Status = 'skipped';
        results(caseIndex).SkipReason = skipReason;
        reportPaths = write_xddx_backend_reports( ...
            results, hostInfo, options, runInfo);
        continue;
    end

    results(caseIndex).Attempted = true;
    set_docker_mode(currentCase.ExecutionMode);
    overrides = overrides_for_case(currentCase);
    try
        timer = tic;
        output = run_default_case( ...
            simulatorFolder, simulatorScript, overrides);
        results(caseIndex).ElapsedSeconds = toc(timer);
        metrics = xddx_backend_norms(output, reference);
        if isempty(reference)
            reference = output;
            referenceCase = currentCase.Name;
            metrics = xddx_backend_norms(output, reference);
        end
        results(caseIndex) = add_success( ...
            results(caseIndex), output, metrics, referenceCase);
    catch exception
        results(caseIndex).Status = 'failed';
        results(caseIndex).ErrorIdentifier = exception.identifier;
        results(caseIndex).ErrorMessage = exception.message;
    end
    reportPaths = write_xddx_backend_reports( ...
        results, hostInfo, options, runInfo);
end
end

function output = run_default_case( ...
    simulatorFolder, simulatorScript, simulatorInputOverrides)
if ~isstruct(simulatorInputOverrides)
    error('xDDx:Tests:InvalidSimulatorOverrides', ...
        'Simulator overrides must be a structure.');
end
oldFolder = cd(simulatorFolder);
folderCleanup = onCleanup(@() cd(oldFolder));
run(simulatorScript);
if ~exist('simulatorData', 'var') ...
        || ~isfield(simulatorData, 'complexPressure')
    error('xDDx:Tests:MissingSimulatorOutput', ...
        'xDDx_simulator did not return simulatorData.complexPressure.');
end
output = simulatorData.complexPressure;
end

function overrides = overrides_for_case(backendCase)
overrides = struct( ...
    'kWaveCalculationFlag', backendCase.Device, ...
    'xDDxCalculationFlag', backendCase.Device, ...
    'cpuArchitecture', 'auto', ...
    'cudaVersion', 'auto', ...
    'useGUI', false, ...
    'shouldPlot', false);
if strcmp(backendCase.Device, 'cpu')
    overrides.cpuArchitecture = backendCase.Variant;
else
    overrides.cudaVersion = backendCase.Variant;
end
end

function set_docker_mode(executionMode)
enabled = char(string(strcmp(executionMode, 'docker')));
setenv('XDDX_USE_DOCKER_ON_WINDOWS', enabled);
setenv('KWAVE_USE_DOCKER_ON_WINDOWS', enabled);
end

function restore_environment(xddxValue, kwaveValue)
setenv('XDDX_USE_DOCKER_ON_WINDOWS', xddxValue);
setenv('KWAVE_USE_DOCKER_ON_WINDOWS', kwaveValue);
end

function result = result_for_case(backendCase)
result = empty_result();
result.Name = backendCase.Name;
result.Device = backendCase.Device;
result.Variant = backendCase.Variant;
result.ExecutionMode = backendCase.ExecutionMode;
end

function result = add_success( ...
    result, output, metrics, referenceCase)
result.Status = 'passed';
result.Success = true;
result.OutputElementCount = numel(output);
result.OutputSize = char(strjoin(string(size(output)), 'x'));
result.OutputL2 = metrics.OutputL2;
result.OutputLInf = metrics.OutputLInf;
result.DifferenceL2 = metrics.DifferenceL2;
result.DifferenceLInf = metrics.DifferenceLInf;
result.RelativeL2 = metrics.RelativeL2;
result.RelativeLInf = metrics.RelativeLInf;
result.ReferenceCase = referenceCase;
end

function value = empty_result()
value = struct( ...
    'Name', '', ...
    'Device', '', ...
    'Variant', '', ...
    'ExecutionMode', '', ...
    'Status', 'pending', ...
    'Attempted', false, ...
    'Success', false, ...
    'SkipReason', '', ...
    'ElapsedSeconds', NaN, ...
    'OutputElementCount', NaN, ...
    'OutputSize', '', ...
    'OutputL2', NaN, ...
    'OutputLInf', NaN, ...
    'DifferenceL2', NaN, ...
    'DifferenceLInf', NaN, ...
    'RelativeL2', NaN, ...
    'RelativeLInf', NaN, ...
    'ReferenceCase', '', ...
    'ErrorIdentifier', '', ...
    'ErrorMessage', '');
end

function runInfo = create_run_info()
timestamp = datetime('now', 'TimeZone', 'UTC');
timestamp.Format = 'yyyyMMdd''T''HHmmss';
runInfo = struct();
runInfo.RunId = char(timestamp);
timestamp.Format = 'yyyy-MM-dd''T''HH:mm:ssXXX';
runInfo.StartedUtc = char(timestamp);
runInfo.CaseDefinition = ...
    'xDDx_simulator defaults; backend explicit; GUI and plots disabled';
runInfo.ReportSchemaVersion = 1;
end
