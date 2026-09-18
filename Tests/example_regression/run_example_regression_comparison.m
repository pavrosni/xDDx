function [summaryTable, comparisonTable, speedTable, results] = run_example_regression_comparison(varargin)
%RUN_EXAMPLE_REGRESSION_COMPARISON Compare example outputs between xDDx versions.
%
% Examples:
%   list_example_regression_scripts
%   run_example_regression_comparison('Scripts', 'fp_sf')
%   run_example_regression_comparison('Scripts', {'quick_start_flat', 'transducer_simulation_sf'})
%   run_example_regression_comparison('Scripts', 'all', 'SimulationDevices', 'both')
%   run_example_regression_comparison('ReferenceProject', 'oldest')
%   run_example_regression_comparison('RunsPerScript', 3, 'StopOnError', false)
%   run_example_regression_comparison('SimulationDevices', 'both')
%   run_example_regression_comparison('NewExecutionMode', 'docker')
%   run_example_regression_comparison('RunMode', 'new-only')
%   run_example_regression_comparison('SaveReports', false)
%   run_example_regression_comparison('SaveMatFile', true)
%
% Outputs:
%   summaryTable    - status and timing for each script in each project
%   comparisonTable - output norms or relative L2 comparison of shared numeric outputs
%   speedTable      - old/new timing ratios per script and rayleigh_simulator core calls
%   results         - raw result structure saved to the report MAT file

options = parse_options(varargin{:});

if ~is_folder(options.NewProjectDir)
    error('xDDx:ExampleRegression:MissingNewProject', ...
        'NewProjectDir does not exist: %s', options.NewProjectDir);
end

runMode = resolve_run_mode(options.RunMode);
if strcmp(runMode, 'compare') && ~is_folder(options.OldProjectDir)
    error('xDDx:ExampleRegression:MissingOldProject', ...
        'Reference project directory does not exist: %s', options.OldProjectDir);
end
scriptList = select_catalog_entries(example_regression_catalog(), ...
    options.Scripts, options.Toolboxes, options.NewProjectDir, ...
    options.OldProjectDir, runMode);
simulationDevices = resolve_simulation_devices(options.SimulationDevices);
newExecutionMode = resolve_new_execution_mode(options.NewExecutionMode);

if isempty(scriptList)
    error('xDDx:ExampleRegression:NoScripts', ...
        'No matching scripts were found for the selected filters.');
end

if ~is_folder(options.OutputDir)
    mkdir(options.OutputDir);
end

fprintf('Comparing %d example script(s).\n', numel(scriptList));
fprintf('Run mode: %s\n', runMode);
fprintf('New project: %s\n', options.NewProjectDir);
if strcmp(runMode, 'compare')
    fprintf('Reference project (%s): %s\n\n', ...
        options.ReferenceProject, options.OldProjectDir);
end
fprintf('New project execution mode: %s\n\n', newExecutionMode);

results = repmat(empty_pair_result(), numel(scriptList)*numel(simulationDevices), 1);

resultIndex = 0;
for iScript = 1:numel(scriptList)
    item = scriptList(iScript);
    for iDevice = 1:numel(simulationDevices)
        resultIndex = resultIndex + 1;
        simulationDevice = simulationDevices{iDevice};
        fprintf('[%d/%d] %s [%s]\n', resultIndex, numel(results), ...
            item.NewRelativePath, simulationDevice);

        results(resultIndex).relativePath = item.NewRelativePath;
        results(resultIndex).runMode = runMode;
        results(resultIndex).simulationDevice = simulationDevice;
        results(resultIndex).newExecutionMode = newExecutionMode;
        results(resultIndex).referenceProject = options.ReferenceProject;
        results(resultIndex).referenceProjectDir = options.OldProjectDir;
        results(resultIndex).new = run_one_side(options.NewProjectDir, ...
            item.NewRelativePath, item.OutputNames, options, ...
            simulationDevice, newExecutionMode);
        if strcmp(runMode, 'compare')
            results(resultIndex).old = run_one_side(options.OldProjectDir, ...
                item.OldRelativePath, item.OutputNames, options, ...
                simulationDevice, 'legacy');
            results(resultIndex).comparison = compare_output_maps( ...
                results(resultIndex).new.outputMap, results(resultIndex).old.outputMap);
        else
            results(resultIndex).comparison = summarize_new_output_map(results(resultIndex).new.outputMap);
        end

        print_pair_summary(results(resultIndex));

        if options.StopOnError && did_result_fail(results(resultIndex))
            error('xDDx:ExampleRegression:ScriptFailed', ...
                'Stopping after failure in %s [%s].', ...
                item.NewRelativePath, simulationDevice);
        end
    end
end

[summaryTable, comparisonTable, speedTable] = build_report_tables(results);

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
matFile = fullfile(options.OutputDir, ['example_regression_' timestamp '.mat']);
summaryCsv = fullfile(options.OutputDir, ['example_regression_summary_' timestamp '.csv']);
comparisonCsv = fullfile(options.OutputDir, ['example_regression_comparison_' timestamp '.csv']);
speedCsv = fullfile(options.OutputDir, ['example_regression_speed_' timestamp '.csv']);

if options.SaveReports
    writetable(summaryTable, summaryCsv);
    writetable(comparisonTable, comparisonCsv);
    writetable(speedTable, speedCsv);

    if options.SaveMatFile
        save(matFile, 'options', 'results', 'summaryTable', 'comparisonTable', 'speedTable', '-v7.3');
    end

    fprintf('\nSaved report files:\n');
    fprintf('  %s\n', summaryCsv);
    fprintf('  %s\n', comparisonCsv);
    fprintf('  %s\n', speedCsv);
    if options.SaveMatFile
        fprintf('  %s\n', matFile);
    end
else
    fprintf('\nSaveReports=false, so no report files were written.\n');
end
end

function options = parse_options(varargin)
repoDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

parser = inputParser;
parser.FunctionName = 'run_example_regression_comparison';
addParameter(parser, 'NewProjectDir', repoDir, @(x) is_text_scalar(x));
addParameter(parser, 'ReferenceProject', 'xddx', @(x) is_text_scalar(x));
addParameter(parser, 'OldProjectDir', '', @(x) is_text_scalar(x));
addParameter(parser, 'Toolboxes', {'holography_toolbox', 'simulation_toolbox'}, @(x) is_text_array(x));
addParameter(parser, 'Scripts', {}, @(x) is_text_array(x));
addParameter(parser, 'OutputDir', fullfile(repoDir, 'Tests', 'example_regression', 'reports'), @(x) is_text_scalar(x));
addParameter(parser, 'RunsPerScript', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'WarmupRuns', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'SimulationDevices', 'script', @(x) is_text_array(x));
addParameter(parser, 'NewExecutionMode', 'auto', @(x) is_text_scalar(x));
addParameter(parser, 'RunMode', 'auto', @(x) is_text_scalar(x));
addParameter(parser, 'StopOnError', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'CloseFigures', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'FigureVisibility', 'off', @(x) any(strcmpi(char(x), {'on', 'off'})));
addParameter(parser, 'SaveReports', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'SaveMatFile', false, @(x) islogical(x) && isscalar(x));
parse(parser, varargin{:});

options = parser.Results;
options.NewProjectDir = char(options.NewProjectDir);
[options.ReferenceProject, options.OldProjectDir] = resolve_reference_project( ...
    options.ReferenceProject, options.OldProjectDir, repoDir);
options.OutputDir = char(options.OutputDir);
options.Toolboxes = normalize_text_array(options.Toolboxes);
options.Scripts = normalize_text_array(options.Scripts);
options.SimulationDevices = normalize_text_array(options.SimulationDevices);
options.NewExecutionMode = char(options.NewExecutionMode);
options.RunMode = char(options.RunMode);
options.RunsPerScript = double(options.RunsPerScript);
options.WarmupRuns = double(options.WarmupRuns);
options.FigureVisibility = char(options.FigureVisibility);
end

function [referenceProject, referenceDir] = resolve_reference_project( ...
        selection, explicitDir, repoDir)
repoParent = fileparts(repoDir);
xddxDir = fullfile(repoParent, 'xDDx');
oldestDir = fullfile(repoParent, 'xDDx_old');
explicitDir = strtrim(char(explicitDir));

if ~isempty(explicitDir)
    referenceDir = explicitDir;
    referenceProject = identify_reference_project(referenceDir, xddxDir, oldestDir);
    return;
end

selection = lower(strtrim(char(selection)));
switch selection
    case {'xddx', 'previous'}
        referenceProject = 'xddx';
        referenceDir = xddxDir;
    case {'oldest', 'xddx_old'}
        referenceProject = 'oldest';
        referenceDir = oldestDir;
    case 'custom'
        error('xDDx:ExampleRegression:MissingCustomReferenceProject', ...
            ['ReferenceProject=''custom'' requires an explicit ' ...
            'OldProjectDir.']);
    otherwise
        error('xDDx:ExampleRegression:InvalidReferenceProject', ...
            ['ReferenceProject must be ''xddx'', ''oldest'', or ''custom''. ' ...
            'Use OldProjectDir to supply a custom directory.']);
end
end

function referenceProject = identify_reference_project(referenceDir, xddxDir, oldestDir)
normalizedReference = normalize_path_for_comparison(referenceDir);
if strcmp(normalizedReference, normalize_path_for_comparison(xddxDir))
    referenceProject = 'xddx';
elseif strcmp(normalizedReference, normalize_path_for_comparison(oldestDir))
    referenceProject = 'oldest';
else
    referenceProject = 'custom';
end
end

function normalizedPath = normalize_path_for_comparison(pathValue)
normalizedPath = strrep(strtrim(char(pathValue)), '/', filesep);
while numel(normalizedPath) > 1 && normalizedPath(end) == filesep
    normalizedPath(end) = [];
end
if ispc
    normalizedPath = lower(normalizedPath);
end
end

function runMode = resolve_run_mode(selection)
selection = lower(strtrim(char(selection)));
allowedModes = {'auto', 'compare', 'new-only'};
if ~any(strcmp(selection, allowedModes))
    error('xDDx:ExampleRegression:InvalidRunMode', ...
        'RunMode must be ''auto'', ''compare'', or ''new-only''.');
end

if ispc
    if strcmp(selection, 'auto')
        runMode = 'compare';
    else
        runMode = selection;
    end
else
    if strcmp(selection, 'compare')
        warning('xDDx:ExampleRegression:NewOnlyForced', ...
            'RunMode=''compare'' was requested, but Linux/macOS cannot run the old project. Using RunMode=''new-only''.');
    end
    runMode = 'new-only';
end
end

function tf = is_text_scalar(value)
tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
    (isstring(value) && isscalar(value));
end

function tf = is_text_array(value)
tf = ischar(value) || isstring(value) || ...
    (iscell(value) && all(cellfun(@is_text_scalar, value)));
end

function values = normalize_text_array(values)
if isempty(values)
    values = {};
elseif ischar(values)
    values = {values};
elseif isstring(values)
    values = cellstr(values(:).');
elseif iscell(values)
    values = cellfun(@char, values(:).', 'UniformOutput', false);
else
    error('xDDx:ExampleRegression:InvalidTextArray', ...
        'Expected text as a character vector, string, or cell array of text.');
end
values = cellfun(@strtrim, values, 'UniformOutput', false);
values = values(~cellfun(@isempty, values));
end

function tf = is_folder(pathValue)
tf = exist(pathValue, 'dir') == 7;
end

function selected = select_catalog_entries(catalog, selection, toolboxNames, ...
        newProjectDir, oldProjectDir, runMode)
catalog = filter_catalog_by_toolbox(catalog, toolboxNames);
selection = normalize_script_selection(selection);
if isempty(selection) || any(strcmpi(selection, 'all'))
    selected = catalog;
else
    selected = catalog(false(size(catalog)));
    for selectionIndex = 1:numel(selection)
        matchIndex = match_catalog_entry(selection{selectionIndex}, catalog);
        selected(end + 1) = catalog(matchIndex); %#ok<AGROW>
    end
    [~, uniqueIndex] = unique({selected.Name}, 'stable');
    selected = selected(sort(uniqueIndex));
end

for itemIndex = 1:numel(selected)
    assert_script_exists(newProjectDir, selected(itemIndex).NewRelativePath, 'new');
    if strcmp(runMode, 'compare')
        assert_script_exists(oldProjectDir, selected(itemIndex).OldRelativePath, 'old');
    end
end
end

function catalog = filter_catalog_by_toolbox(catalog, toolboxNames)
if isempty(toolboxNames) || any(strcmpi(toolboxNames, 'all'))
    return;
end

requested = lower(toolboxNames);
allowed = {'holography_toolbox', 'simulation_toolbox'};
invalid = setdiff(requested, allowed);
if ~isempty(invalid)
    error('xDDx:ExampleRegression:InvalidToolbox', ...
        'Toolboxes must contain holography_toolbox, simulation_toolbox, or all. Invalid: %s', ...
        strjoin(invalid, ', '));
end

paths = lower({catalog.NewRelativePath});
keep = false(size(catalog));
for toolboxIndex = 1:numel(requested)
    prefix = [requested{toolboxIndex} '/'];
    keep = keep | cellfun(@(pathValue) strncmp(pathValue, prefix, numel(prefix)), paths);
end
catalog = catalog(keep);
end

function selection = normalize_script_selection(selection)
selection = normalize_text_array(selection);
selection = cellfun(@(value) strrep(value, '\', '/'), selection, ...
    'UniformOutput', false);
end

function matchIndex = match_catalog_entry(selection, catalog)
selection = remove_m_extension(strtrim(selection));
names = {catalog.Name};
newPaths = cellfun(@remove_m_extension, ...
    {catalog.NewRelativePath}, 'UniformOutput', false);
oldPaths = cellfun(@remove_m_extension, ...
    {catalog.OldRelativePath}, 'UniformOutput', false);
baseNames = cell(size(newPaths));
for pathIndex = 1:numel(newPaths)
    [~, baseNames{pathIndex}] = fileparts(newPaths{pathIndex});
end

isMatch = strcmpi(names, selection) | strcmpi(newPaths, selection) | ...
    strcmpi(oldPaths, selection) | strcmpi(baseNames, selection);
matchIndex = find(isMatch);
if isempty(matchIndex)
    error('xDDx:ExampleRegression:UnknownScript', ...
        ['Selected script was not found: %s. Run ' ...
        'list_example_regression_scripts to see valid selections.'], selection);
elseif numel(matchIndex) > 1
    error('xDDx:ExampleRegression:AmbiguousScript', ...
        'Script selection is ambiguous; use a relative path: %s', selection);
end
end

function assert_script_exists(projectDir, relativePath, projectLabel)
absolutePath = fullfile(projectDir, relativePath);
if exist(absolutePath, 'file') ~= 2
    error('xDDx:ExampleRegression:MissingScript', ...
        'Mapped %s script does not exist: %s', projectLabel, absolutePath);
end
end

function simulationDevices = resolve_simulation_devices(selection)
selection = lower_cellstr(selection);

if isempty(selection) || any(strcmp(selection, 'script'))
    simulationDevices = {'script'};
    return;
end

if any(strcmp(selection, 'both'))
    simulationDevices = {'cpu'};
    if is_cuda_available()
        simulationDevices = {'cuda', 'cpu'};
    else
        warning('xDDx:ExampleRegression:CudaUnavailable', ...
            'SimulationDevices=''both'' was requested, but CUDA was not detected. Running CPU only.');
    end
else
    allowedDevices = {'cuda', 'cpu'};
    invalidDevices = setdiff(selection, allowedDevices);
    if ~isempty(invalidDevices)
        error('xDDx:ExampleRegression:InvalidSimulationDevice', ...
            'SimulationDevices must be ''script'', ''both'', ''cuda'', ''cpu'', or a cell array of ''cuda''/''cpu''. Invalid: %s', ...
            strjoin(invalidDevices, ', '));
    end

    simulationDevices = unique(selection, 'stable');
    if any(strcmp(simulationDevices, 'cuda')) && ~is_cuda_available()
        warning('xDDx:ExampleRegression:CudaUnavailable', ...
            'CUDA was requested but not detected. CUDA runs may fail.');
    end
end
end

function executionMode = resolve_new_execution_mode(selection)
selection = lower(strtrim(char(selection)));

if ispc
    if strcmp(selection, 'auto')
        executionMode = 'native';
    elseif any(strcmp(selection, {'native', 'docker'}))
        executionMode = selection;
    else
        error('xDDx:ExampleRegression:InvalidNewExecutionMode', ...
            'NewExecutionMode must be ''auto'', ''native'', or ''docker''.');
    end
else
    if strcmp(selection, 'native')
        warning('xDDx:ExampleRegression:DockerForced', ...
            'NewExecutionMode=''native'' was requested, but Linux/macOS tests must use Docker. Using Docker.');
    elseif ~any(strcmp(selection, {'auto', 'docker'}))
        error('xDDx:ExampleRegression:InvalidNewExecutionMode', ...
            'NewExecutionMode must be ''auto'', ''native'', or ''docker''.');
    end
    executionMode = 'docker';
end
end

function values = lower_cellstr(values)
values = cellfun(@lower, values, 'UniformOutput', false);
end

function tf = is_cuda_available()
tf = is_cuda_available_from_gpu_device_count();
if tf
    return;
end

tf = is_cuda_available_from_gpu_device();
if tf
    return;
end

tf = is_cuda_available_from_nvidia_smi();
end

function tf = is_cuda_available_from_gpu_device_count()
tf = false;
if exist('gpuDeviceCount', 'file') == 2
    try
        tf = gpuDeviceCount('available') > 0;
    catch
        tf = false;
    end
end
end

function tf = is_cuda_available_from_gpu_device()
tf = false;
if exist('gpuDevice', 'file') == 2
    try
        gpuDevice;
        tf = true;
    catch
        tf = false;
    end
end
end

function tf = is_cuda_available_from_nvidia_smi()
try
    [status, output] = system('nvidia-smi -L');
    tf = status == 0 && contains(output, 'GPU');
catch
    tf = false;
end
end

function value = remove_m_extension(value)
if numel(value) >= 2 && strcmpi(value((end - 1):end), '.m')
    value = value(1:(end - 2));
end
end

function pairResult = empty_pair_result()
pairResult = struct();
pairResult.relativePath = '';
pairResult.runMode = '';
pairResult.simulationDevice = '';
pairResult.newExecutionMode = '';
pairResult.referenceProject = '';
pairResult.referenceProjectDir = '';
pairResult.new = empty_side_result();
pairResult.old = empty_side_result();
pairResult.comparison = struct([]);
end

function tf = did_result_fail(result)
tf = ~result.new.success;
if strcmp(result.runMode, 'compare')
    tf = tf || ~result.old.success;
end
end

function sideResult = empty_side_result()
sideResult = struct();
sideResult.success = false;
sideResult.elapsedSeconds = NaN;
sideResult.elapsedSecondsAllRuns = [];
sideResult.coreTotalSeconds = NaN;
sideResult.coreTotalSecondsAllRuns = [];
sideResult.coreCalls = NaN;
sideResult.coreCallsAllRuns = [];
sideResult.coreMeanSecondsPerCall = NaN;
sideResult.coreMeanSecondsPerCallAllRuns = [];
sideResult.executionMode = '';
sideResult.errorIdentifier = '';
sideResult.errorMessage = '';
sideResult.expectedOutputNames = {};
sideResult.missingOutputNames = {};
sideResult.outputMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

function result = run_one_side(projectDir, relativePath, outputNames, ...
        options, simulationDevice, executionMode)
result = empty_side_result();
result.executionMode = executionMode;
scriptPath = fullfile(projectDir, relativePath);
scriptFolder = fileparts(scriptPath);
[scriptRunCommand, tempScriptPath, generatedFiles] = ...
    build_script_run_command(scriptPath, simulationDevice);

oldPath = path;
oldFolder = pwd;
oldDefaultVisibility = get(groot, 'DefaultFigureVisible');
oldDockerEnv = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
cleanup = onCleanup(@() restore_session(oldPath, oldFolder, ...
    oldDefaultVisibility, oldDockerEnv, tempScriptPath, ...
    generatedFiles, options.CloseFigures));

path(pathdef);
clear_functions();
cd(scriptFolder);
set(groot, 'DefaultFigureVisible', options.FigureVisibility);
apply_execution_mode(executionMode);

try
    for iWarmup = 1:options.WarmupRuns
        close_figures_if_requested(options.CloseFigures);
        evalin('base', 'clear variables;');
        evalin('base', scriptRunCommand);
    end

    elapsed = zeros(1, options.RunsPerScript);
    coreTotalSeconds = NaN(1, options.RunsPerScript);
    coreCalls = NaN(1, options.RunsPerScript);
    coreMeanSecondsPerCall = NaN(1, options.RunsPerScript);
    outputMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for iRun = 1:options.RunsPerScript
        close_figures_if_requested(options.CloseFigures);
        evalin('base', 'clear variables;');
        profile clear;
        profile on;
        timerHandle = tic;
        evalin('base', scriptRunCommand);
        elapsed(iRun) = toc(timerHandle);
        profile off;
        coreInfo = extract_rayleigh_profile_info(profile('info'));
        coreTotalSeconds(iRun) = coreInfo.totalSeconds;
        coreCalls(iRun) = coreInfo.calls;
        coreMeanSecondsPerCall(iRun) = coreInfo.meanSecondsPerCall;
        if iRun == options.RunsPerScript
            [outputMap, expectedOutputNames, missingOutputNames] = ...
                capture_numeric_outputs(outputNames);
        end
    end

    result.success = true;
    result.elapsedSecondsAllRuns = elapsed;
    result.elapsedSeconds = median(elapsed);
    result.coreTotalSecondsAllRuns = coreTotalSeconds;
    result.coreTotalSeconds = median_omitnan(coreTotalSeconds);
    result.coreCallsAllRuns = coreCalls;
    result.coreCalls = median_omitnan(coreCalls);
    result.coreMeanSecondsPerCallAllRuns = coreMeanSecondsPerCall;
    result.coreMeanSecondsPerCall = median_omitnan(coreMeanSecondsPerCall);
    result.expectedOutputNames = expectedOutputNames;
    result.missingOutputNames = missingOutputNames;
    result.outputMap = outputMap;
catch ME
    profile off;
    result.success = false;
    result.errorIdentifier = ME.identifier;
    result.errorMessage = ME.message;
end
end

function value = median_omitnan(values)
values = values(~isnan(values));
if isempty(values)
    value = NaN;
else
    value = median(values);
end
end

function apply_execution_mode(executionMode)
switch executionMode
    case 'docker'
        setenv('XDDX_USE_DOCKER_ON_WINDOWS', '1');
    case 'native'
        setenv('XDDX_USE_DOCKER_ON_WINDOWS', '0');
end
end

function [scriptRunCommand, tempScriptPath, generatedFiles] = ...
        build_script_run_command(scriptPath, simulationDevice)
generatedFiles = {fullfile(tempdir, 'xddx_example_regression_output.mat')};
scriptBytes = read_binary_file(scriptPath);
assignmentMatch = '';
if ~strcmp(simulationDevice, 'script')
    replacement = sprintf('simulationDevice = ''%s'';', simulationDevice);
    assignmentPattern = 'simulationDevice\s*=\s*[''"][^''"]+[''"]\s*;';
    [scriptBytes, assignmentMatch] = replace_assignment_in_bytes( ...
        scriptBytes, assignmentPattern, replacement);
end

outputReplacement = sprintf('outFilename = ''%s'';', ...
    escape_matlab_string(generatedFiles{1}));
outputPattern = 'outFilename\s*=\s*[''"][^''"]+[''"]\s*;';
[scriptBytes, ~] = replace_assignment_in_bytes( ...
    scriptBytes, outputPattern, outputReplacement);
scriptBytes = prepend_script_folder_change(scriptBytes, fileparts(scriptPath));

scriptText = char(scriptBytes(:).');
if ~strcmp(simulationDevice, 'script') && isempty(assignmentMatch) ...
        && contains(scriptText, 'rayleigh_simulator')
    warning('xDDx:ExampleRegression:NoSimulationDeviceAssignment', ...
        'No simulationDevice assignment found in %s. Running without device override.', scriptPath);
end

tempScriptPath = make_temp_script_path();
write_binary_file(tempScriptPath, scriptBytes);
scriptRunCommand = sprintf('run(''%s'');', escape_matlab_string(tempScriptPath));
end

function scriptBytes = prepend_script_folder_change(scriptBytes, scriptFolder)
scriptPrelude = sprintf('cd(''%s'');\n', escape_matlab_string(scriptFolder));
scriptBytes = [uint8(scriptPrelude(:)); scriptBytes];
end

function [scriptBytes, assignmentMatch] = replace_assignment_in_bytes( ...
        scriptBytes, assignmentPattern, replacement)
scriptTextForAsciiSearch = char(scriptBytes(:).');
[startIndex, endIndex, assignmentMatch] = regexp(scriptTextForAsciiSearch, assignmentPattern, 'start', 'end', 'match', 'once');

if isempty(assignmentMatch)
    return;
end

scriptBytes = [scriptBytes(1:(startIndex - 1)); uint8(replacement(:)); scriptBytes((endIndex + 1):end)];
end

function bytes = read_binary_file(filePath)
fid = fopen(filePath, 'r');
if fid < 0
    error('xDDx:ExampleRegression:ScriptReadFailed', ...
        'Could not read script: %s', filePath);
end
cleanupFile = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
end

function tempScriptPath = make_temp_script_path()
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
baseName = sprintf('xddx_regression_tmp_%s_%06d.m', ...
    timestamp, randi(999999));
tempScriptPath = fullfile(tempdir, baseName);
end

function write_binary_file(filePath, bytes)
fid = fopen(filePath, 'w');
if fid < 0
    error('xDDx:ExampleRegression:TempScriptWriteFailed', ...
        'Could not create temporary script: %s', filePath);
end
cleanupFile = onCleanup(@() fclose(fid));
fwrite(fid, bytes, 'uint8');
end

function coreInfo = extract_rayleigh_profile_info(profileInfo)
coreInfo = struct();
coreInfo.totalSeconds = NaN;
coreInfo.calls = 0;
coreInfo.meanSecondsPerCall = NaN;

if ~isfield(profileInfo, 'FunctionTable') || isempty(profileInfo.FunctionTable)
    return;
end

functionTable = profileInfo.FunctionTable;
isRayleighSimulator = false(1, numel(functionTable));
for iFunction = 1:numel(functionTable)
    functionName = functionTable(iFunction).FunctionName;
    if contains(functionName, '>')
        continue;
    end
    [~, baseName] = fileparts(char(functionName));
    isRayleighSimulator(iFunction) = strcmp(baseName, 'rayleigh_simulator');
end

if ~any(isRayleighSimulator)
    return;
end

coreRows = functionTable(isRayleighSimulator);
coreInfo.totalSeconds = sum([coreRows.TotalTime]);
coreInfo.calls = sum([coreRows.NumCalls]);
if coreInfo.calls > 0
    coreInfo.meanSecondsPerCall = coreInfo.totalSeconds / coreInfo.calls;
end
end

function text = escape_matlab_string(text)
text = char(text);
text = strrep(text, '''', '''''');
end

function restore_session(oldPath, oldFolder, oldDefaultVisibility, ...
        oldDockerEnv, tempScriptPath, generatedFiles, closeFigures)
profile off;
close_figures_if_requested(closeFigures);
evalin('base', 'clear variables;');
delete_temp_script(tempScriptPath);
delete_generated_files(generatedFiles);
path(oldPath);
cd(oldFolder);
set(groot, 'DefaultFigureVisible', oldDefaultVisibility);
setenv('XDDX_USE_DOCKER_ON_WINDOWS', oldDockerEnv);
clear_functions();
end

function delete_generated_files(generatedFiles)
for fileIndex = 1:numel(generatedFiles)
    if exist(generatedFiles{fileIndex}, 'file') == 2
        delete(generatedFiles{fileIndex});
    end
end
end

function delete_temp_script(tempScriptPath)
if ~isempty(tempScriptPath) && exist(tempScriptPath, 'file') == 2
    delete(tempScriptPath);
end
end

function close_figures_if_requested(closeFigures)
if closeFigures
    close all force;
end
end

function clear_functions()
% Script runs use isolated paths/workspaces; no function cache reset needed.
end

function [outputMap, expectedOutputNames, missingOutputNames] = ...
        capture_numeric_outputs(expectedOutputNames)
outputMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
missingOutputNames = {};

for iOutput = 1:numel(expectedOutputNames)
    outputName = expectedOutputNames{iOutput};
    [exists, value] = try_get_output_value(outputName);
    if exists
        if isnumeric(value) || islogical(value)
            if ~isempty(value) && all(isfinite(double(value(:))) | isnan(double(value(:))))
                outputMap(char(outputName)) = value;
            end
        end
    else
        missingOutputNames(end + 1) = outputName; %#ok<AGROW>
    end
end
end

function [exists, value] = try_get_output_value(outputName)
value = [];
try
    value = evalin('base', char(outputName));
    exists = true;
catch
    exists = false;
end
end

function outputMap = flatten_numeric_value(value, prefix)
outputMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

if isnumeric(value) || islogical(value)
    if ~isempty(value) && all(isfinite(double(value(:))) | isnan(double(value(:))))
        outputMap(prefix) = value;
    end
elseif isstruct(value) && isscalar(value)
    fieldNames = fieldnames(value);
    for iField = 1:numel(fieldNames)
        childMap = flatten_numeric_value(value.(fieldNames{iField}), ...
            sprintf('%s.%s', prefix, fieldNames{iField}));
        outputMap = merge_maps(outputMap, childMap);
    end
elseif iscell(value)
    for iCell = 1:numel(value)
        childMap = flatten_numeric_value(value{iCell}, sprintf('%s{%d}', prefix, iCell));
        outputMap = merge_maps(outputMap, childMap);
    end
end
end

function target = merge_maps(target, source)
keys = source.keys;
for iKey = 1:numel(keys)
    target(keys{iKey}) = source(keys{iKey});
end
end

function comparisons = compare_output_maps(newMap, oldMap)
newKeys = newMap.keys;
oldKeys = oldMap.keys;
sharedKeys = intersect(newKeys, oldKeys, 'stable');

comparisons = repmat(empty_comparison(), numel(sharedKeys), 1);

for iKey = 1:numel(sharedKeys)
    key = char(sharedKeys(iKey));
    newValue = newMap(key);
    oldValue = oldMap(key);

    comparisons(iKey).outputName = key;
    comparisons(iKey).newSize = mat2str(size(newValue));
    comparisons(iKey).oldSize = mat2str(size(oldValue));
    comparisons(iKey).sameSize = isequal(size(newValue), size(oldValue));

    if comparisons(iKey).sameSize
        newVector = double(newValue(:));
        oldVector = double(oldValue(:));
        difference = newVector - oldVector;
        denominator = max(norm(oldVector, 2), eps);

        comparisons(iKey).relativeL2 = norm(difference, 2) / denominator;
        comparisons(iKey).absoluteL2 = norm(difference, 2);
        comparisons(iKey).oldL2 = norm(oldVector, 2);
        comparisons(iKey).maxAbsDifference = max_omitnan(abs(difference));
    end
end
end

function comparisons = summarize_new_output_map(newMap)
newKeys = newMap.keys;
comparisons = repmat(empty_comparison(), numel(newKeys), 1);

for iKey = 1:numel(newKeys)
    key = char(newKeys(iKey));
    newValue = newMap(key);
    comparisons(iKey).outputName = key;
    comparisons(iKey).newSize = mat2str(size(newValue));
    comparisons(iKey).newL2 = norm(double(newValue(:)), 2);
end
end

function value = max_omitnan(values)
values = values(~isnan(values));
if isempty(values)
    value = NaN;
else
    value = max(values);
end
end

function comparison = empty_comparison()
comparison = struct();
comparison.outputName = '';
comparison.newSize = '';
comparison.newL2 = NaN;
comparison.oldSize = '';
comparison.sameSize = false;
comparison.relativeL2 = NaN;
comparison.absoluteL2 = NaN;
comparison.oldL2 = NaN;
comparison.maxAbsDifference = NaN;
end

function [summaryTable, comparisonTable, speedTable] = build_report_tables(results)
summaryRows = cell(0, 18);
row = 0;
for iResult = 1:numel(results)
    scriptLabel = short_script_name(results(iResult).relativePath);
    deviceLabel = char(results(iResult).simulationDevice);
    row = row + 1;
    summaryRows(row, :) = make_summary_row(scriptLabel, deviceLabel, 'new', results(iResult).new);
    if strcmp(results(iResult).runMode, 'compare')
        row = row + 1;
        summaryRows(row, :) = make_summary_row(scriptLabel, deviceLabel, ...
            results(iResult).referenceProject, results(iResult).old);
    end
end

summaryTable = cell2table(summaryRows, 'VariableNames', { ...
    'Script', 'Device', 'Project', 'ExecutionMode', 'Success', 'ElapsedSecondsMedian', 'ElapsedSecondsAllRuns', ...
    'CoreTotalSecondsMedian', 'CoreTotalSecondsAllRuns', 'CoreCallsMedian', ...
    'CoreCallsAllRuns', 'CoreMeanSecondsPerCall', 'CoreMeanSecondsPerCallAllRuns', ...
    'OutputCount', 'ExpectedOutputNames', 'MissingOutputNames', 'ErrorIdentifier', 'ErrorMessage'});

if all(strcmp({results.runMode}, 'new-only'))
    comparisonTable = build_new_only_output_table(results);
else
    comparisonTable = build_comparison_table(results);
end

speedRows = cell(numel(results), 18);
for iResult = 1:numel(results)
    newTime = results(iResult).new.elapsedSeconds;
    newCoreTotal = results(iResult).new.coreTotalSeconds;
    newCoreCalls = results(iResult).new.coreCalls;
    newCoreMean = results(iResult).new.coreMeanSecondsPerCall;
    oldTime = NaN;
    oldCoreTotal = NaN;
    oldCoreCalls = NaN;
    oldCoreMean = NaN;
    bothSucceeded = results(iResult).new.success;
    if strcmp(results(iResult).runMode, 'compare')
        oldTime = results(iResult).old.elapsedSeconds;
        oldCoreTotal = results(iResult).old.coreTotalSeconds;
        oldCoreCalls = results(iResult).old.coreCalls;
        oldCoreMean = results(iResult).old.coreMeanSecondsPerCall;
        bothSucceeded = results(iResult).new.success && results(iResult).old.success;
    end
    speedRows(iResult, :) = {short_script_name(results(iResult).relativePath), ...
        char(results(iResult).simulationDevice), ...
        char(results(iResult).runMode), ...
        char(results(iResult).newExecutionMode), ...
        bothSucceeded, ...
        newTime, oldTime, safe_ratio(oldTime, newTime), safe_ratio(newTime, oldTime), ...
        safe_relative_difference(newTime, oldTime), ...
        newCoreTotal, oldCoreTotal, safe_ratio(oldCoreTotal, newCoreTotal), ...
        newCoreCalls, oldCoreCalls, newCoreMean, oldCoreMean, ...
        safe_relative_difference(newCoreMean, oldCoreMean)};
end

speedTable = cell2table(speedRows, 'VariableNames', { ...
    'Script', 'Device', 'RunMode', 'NewExecutionMode', 'BothSucceeded', 'NewElapsedSecondsMedian', ...
    'OldElapsedSecondsMedian', 'OldOverNewTimeRatio', 'NewOverOldTimeRatio', ...
    'ElapsedRelativeDifferenceNewVsOld', ...
    'NewCoreTotalSecondsMedian', 'OldCoreTotalSecondsMedian', ...
    'OldOverNewCoreTotalTimeRatio', 'NewCoreCallsMedian', 'OldCoreCallsMedian', ...
    'NewCoreMeanSecondsPerCall', 'OldCoreMeanSecondsPerCall', ...
    'CoreMeanRelativeDifferenceNewVsOld'});
end

function outputTable = build_new_only_output_table(results)
outputRows = {};
for iResult = 1:numel(results)
    scriptLabel = short_script_name(results(iResult).relativePath);
    deviceLabel = char(results(iResult).simulationDevice);
    runModeLabel = char(results(iResult).runMode);
    executionModeLabel = char(results(iResult).newExecutionMode);
    outputs = results(iResult).comparison;
    if isempty(outputs)
        outputRows(end + 1, :) = {scriptLabel, deviceLabel, runModeLabel, executionModeLabel, '', '', NaN}; %#ok<AGROW>
        continue;
    end

    for iOutput = 1:numel(outputs)
        outputRows(end + 1, :) = {scriptLabel, deviceLabel, runModeLabel, executionModeLabel, ...
            outputs(iOutput).outputName, outputs(iOutput).newSize, outputs(iOutput).newL2}; %#ok<AGROW>
    end
end

outputTable = cell2table(outputRows, 'VariableNames', { ...
    'Script', 'Device', 'RunMode', 'NewExecutionMode', 'OutputName', 'NewSize', 'NewL2'});
end

function comparisonTable = build_comparison_table(results)
comparisonRows = {};
for iResult = 1:numel(results)
    scriptLabel = short_script_name(results(iResult).relativePath);
    deviceLabel = char(results(iResult).simulationDevice);
    executionModeLabel = char(results(iResult).newExecutionMode);
    runModeLabel = char(results(iResult).runMode);
    comparisons = results(iResult).comparison;
    if isempty(comparisons)
        comparisonRows(end + 1, :) = {scriptLabel, deviceLabel, runModeLabel, executionModeLabel, '', '', '', ...
            false, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        continue;
    end

    for iComparison = 1:numel(comparisons)
        c = comparisons(iComparison);
        comparisonRows(end + 1, :) = {scriptLabel, deviceLabel, runModeLabel, executionModeLabel, c.outputName, ...
            c.newSize, c.oldSize, c.sameSize, c.relativeL2, c.absoluteL2, ...
            c.oldL2, c.maxAbsDifference}; %#ok<AGROW>
    end
end

comparisonTable = cell2table(comparisonRows, 'VariableNames', { ...
    'Script', 'Device', 'RunMode', 'NewExecutionMode', 'OutputName', 'NewSize', 'OldSize', 'SameSize', ...
    'RelativeL2', 'AbsoluteL2', 'OldL2', 'MaxAbsDifference'});
end

function ratio = safe_ratio(numerator, denominator)
if isnan(numerator) || isnan(denominator) || denominator == 0
    ratio = NaN;
else
    ratio = numerator / denominator;
end
end

function relativeDifference = safe_relative_difference(newValue, oldValue)
if isnan(newValue) || isnan(oldValue) || oldValue == 0
    relativeDifference = NaN;
else
    relativeDifference = (newValue - oldValue) / oldValue;
end
end

function scriptLabel = short_script_name(relativePath)
[~, name, ext] = fileparts(char(relativePath));
scriptLabel = [name ext];
end

function row = make_summary_row(relativePath, deviceLabel, projectLabel, sideResult)
row = {relativePath, deviceLabel, projectLabel, char(sideResult.executionMode), sideResult.success, sideResult.elapsedSeconds, ...
    mat2str(sideResult.elapsedSecondsAllRuns), ...
    sideResult.coreTotalSeconds, mat2str(sideResult.coreTotalSecondsAllRuns), ...
    sideResult.coreCalls, mat2str(sideResult.coreCallsAllRuns), ...
    sideResult.coreMeanSecondsPerCall, mat2str(sideResult.coreMeanSecondsPerCallAllRuns), ...
    sideResult.outputMap.Count, ...
    strjoin(sideResult.expectedOutputNames, '; '), strjoin(sideResult.missingOutputNames, '; '), ...
    char(sideResult.errorIdentifier), char(sideResult.errorMessage)};
end

function print_pair_summary(result)
newStatus = 'failed';
oldStatus = 'failed';
if result.new.success
    newStatus = sprintf('ok %.3f s', result.new.elapsedSeconds);
end

fprintf('  new: %s\n', newStatus);
if strcmp(result.runMode, 'compare')
    if result.old.success
        oldStatus = sprintf('ok %.3f s', result.old.elapsedSeconds);
    end
    fprintf('  old: %s\n', oldStatus);
end

comparisons = result.comparison;
sameSize = [comparisons.sameSize];
relativeL2 = [comparisons.relativeL2];
if ~strcmp(result.runMode, 'compare')
    fprintf('  curated outputs captured: %d\n\n', numel(comparisons));
elseif any(sameSize)
    fprintf('  curated outputs: %d, same-size comparable: %d, max relative L2: %.3g\n\n', ...
        numel(comparisons), nnz(sameSize), max_omitnan(relativeL2(sameSize)));
else
    fprintf('  curated outputs: %d, same-size comparable: 0\n\n', numel(comparisons));
end
end
