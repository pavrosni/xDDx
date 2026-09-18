function [simulatorInputs, autoDeviceSelection] = resolve_auto_simulation_devices(simulatorInputs, varargin)
%RESOLVE_AUTO_SIMULATION_DEVICES Resolve and cache simulator device flags.
%
% Use 'auto' for kWaveCalculationFlag and/or xDDxCalculationFlag to select
% CUDA on CUDA-capable Windows/Linux machines and CPU otherwise. macOS is
% always resolved to CPU. The resolved choice is cached in MATLAB settings
% after a successful simulation so repeated runs on the same PC avoid the
% CUDA detection checks.

options = parse_options(varargin{:});

if options.SaveSelection
    save_auto_device_cache(options.AutoDeviceSelection);
    return;
end

if options.ClearCache
    clear_auto_device_cache();
    return;
end

kWaveCalculationFlag = normalize_calculation_flag(simulatorInputs.kWaveCalculationFlag, 'kWaveCalculationFlag');
xDDxCalculationFlag = normalize_calculation_flag(simulatorInputs.xDDxCalculationFlag, 'xDDxCalculationFlag');

autoDeviceSelection = create_empty_selection(kWaveCalculationFlag, xDDxCalculationFlag);
usesAuto = strcmp(kWaveCalculationFlag, 'auto') || strcmp(xDDxCalculationFlag, 'auto');
autoDeviceSelection.usesAuto = usesAuto;

if ~usesAuto
    simulatorInputs.kWaveCalculationFlag = kWaveCalculationFlag;
    simulatorInputs.xDDxCalculationFlag = xDDxCalculationFlag;
    return;
end

cachedSelection = [];
if ~options.ForceDetect
    cachedSelection = load_auto_device_cache();
end

if ~isempty(cachedSelection)
    selectedDevice = cachedSelection.selectedDevice;
    autoDeviceSelection.usedCache = true;
    autoDeviceSelection.detectionSummary = cachedSelection.detectionSummary;
else
    warning('resolve_auto_simulation_devices:AutoDetectionMayTakeTime', ...
        ['Automatic simulation device detection is running without a cached selection. ' ...
        'This may take several minutes, especially with older MATLAB releases.\n' ...
        'If detection appears stuck, see ' ...
        'https://github.com/pavrosni/xDDx/blob/main/docs/BACKEND_SELECTION.md\n' ...
        'To bypass automatic detection on Windows/Linux Intel or AMD x86-64 computers, ' ...
        'edit xDDx_simulator.m and rerun with:\n' ...
        'simulatorInputs.kWaveCalculationFlag = ''cpu'';\n' ...
        'simulatorInputs.xDDxCalculationFlag = ''cpu'';\n' ...
        'simulatorInputs.cpuArchitecture = ''sse2'';']);
    [selectedDevice, detectionSummary] = detect_auto_device();
    autoDeviceSelection.usedCache = false;
    autoDeviceSelection.detectionSummary = detectionSummary;
end

if strcmp(kWaveCalculationFlag, 'auto')
    simulatorInputs.kWaveCalculationFlag = selectedDevice;
else
    simulatorInputs.kWaveCalculationFlag = kWaveCalculationFlag;
end

if strcmp(xDDxCalculationFlag, 'auto')
    simulatorInputs.xDDxCalculationFlag = selectedDevice;
else
    simulatorInputs.xDDxCalculationFlag = xDDxCalculationFlag;
end

autoDeviceSelection.selectedDevice = selectedDevice;
autoDeviceSelection.kWaveCalculationFlag = simulatorInputs.kWaveCalculationFlag;
autoDeviceSelection.xDDxCalculationFlag = simulatorInputs.xDDxCalculationFlag;

fprintf('Auto simulation device: k-Wave = %s, xDDx = %s (%s)\n', ...
    simulatorInputs.kWaveCalculationFlag, simulatorInputs.xDDxCalculationFlag, ...
    autoDeviceSelection.detectionSummary);
end

function options = parse_options(varargin)
options = struct();
options.ForceDetect = false;
options.SaveSelection = false;
options.ClearCache = false;
options.AutoDeviceSelection = [];

if mod(numel(varargin), 2) ~= 0
    error('resolve_auto_simulation_devices:InvalidOptions', ...
        'Options must be provided as name-value pairs.');
end

for idx = 1:2:numel(varargin)
    optionName = char(varargin{idx});
    optionValue = varargin{idx + 1};
    switch lower(optionName)
        case 'forcedetect'
            options.ForceDetect = logical(optionValue);
        case 'saveselection'
            options.SaveSelection = logical(optionValue);
        case 'clearcache'
            options.ClearCache = logical(optionValue);
        case 'autodeviceselection'
            options.AutoDeviceSelection = optionValue;
        otherwise
            error('resolve_auto_simulation_devices:InvalidOption', ...
                'Unknown option "%s".', optionName);
    end
end
end

function calculationFlag = normalize_calculation_flag(calculationFlag, variableName)
if isa(calculationFlag, 'string')
    calculationFlag = char(calculationFlag);
end

if ~ischar(calculationFlag)
    error('resolve_auto_simulation_devices:InvalidFlag', ...
        '%s must be ''auto'', ''cpu'', or ''cuda''.', variableName);
end

calculationFlag = lower(strtrim(calculationFlag));
if ~any(strcmp(calculationFlag, {'auto', 'cpu', 'cuda'}))
    error('resolve_auto_simulation_devices:InvalidFlag', ...
        '%s must be ''auto'', ''cpu'', or ''cuda''.', variableName);
end
end

function selection = create_empty_selection(kWaveCalculationFlag, xDDxCalculationFlag)
selection = struct();
selection.schemaVersion = 1;
selection.usesAuto = false;
selection.usedCache = false;
selection.selectedDevice = '';
selection.kWaveCalculationFlag = kWaveCalculationFlag;
selection.xDDxCalculationFlag = xDDxCalculationFlag;
selection.originalKWaveCalculationFlag = kWaveCalculationFlag;
selection.originalXDDxCalculationFlag = xDDxCalculationFlag;
selection.detectionSummary = '';
selection.machineId = get_machine_id();
selection.platform = computer();
selection.matlabRelease = version('-release');
selection.savedAt = '';
end

function [selectedDevice, detectionSummary] = detect_auto_device()
if ismac
    selectedDevice = 'cpu';
    detectionSummary = 'macOS uses CPU mode';
    return;
end

[hasMatlabGpu, matlabGpuSummary] = detect_matlab_cuda_gpu();
if hasMatlabGpu
    selectedDevice = 'cuda';
    detectionSummary = matlabGpuSummary;
    return;
end

[hasNvidiaSmiGpu, nvidiaSmiSummary] = detect_nvidia_smi_gpu();
if hasNvidiaSmiGpu
    selectedDevice = 'cuda';
    detectionSummary = nvidiaSmiSummary;
    return;
end

selectedDevice = 'cpu';
detectionSummary = ['no CUDA-capable GPU detected; ' matlabGpuSummary '; ' nvidiaSmiSummary];
end

function [hasGpu, summary] = detect_matlab_cuda_gpu()
hasGpu = false;
summary = 'MATLAB GPU check unavailable';

if exist('gpuDeviceCount', 'file') ~= 2
    return;
end

try
    gpuCount = gpuDeviceCount('available');
catch
    try
        gpuCount = gpuDeviceCount();
    catch ME
        summary = ['MATLAB GPU check failed: ' ME.message];
        return;
    end
end

if gpuCount < 1
    summary = 'MATLAB found no available GPU devices';
    return;
end

try
    gpuInfo = gpuDevice();
    if isprop(gpuInfo, 'Name')
        summary = ['MATLAB GPU detected: ' char(gpuInfo.Name)];
    else
        summary = 'MATLAB GPU detected';
    end
catch
    summary = sprintf('MATLAB GPU detected (%d device(s))', gpuCount);
end
hasGpu = true;
end

function [hasGpu, summary] = detect_nvidia_smi_gpu()
hasGpu = false;
summary = 'nvidia-smi check unavailable or no NVIDIA GPU detected';

cmdout = get_nvidia_smi_output('-L');
if isempty(cmdout)
    return;
end

cmdoutLower = lower(cmdout);
if contains(cmdoutLower, 'gpu') || contains(cmdoutLower, 'nvidia')
    hasGpu = true;
    firstLine = get_first_line(strtrim(cmdout));
    summary = ['nvidia-smi detected CUDA-capable GPU: ' firstLine];
end
end

function output = get_nvidia_smi_output(commandArguments)
output = '';
commands = get_nvidia_smi_commands();

for iCommand = 1:numel(commands)
    [status, result] = system([commands{iCommand} ' ' commandArguments ' ' get_shell_error_redirect()]);
    if status == 0 && ~isempty(strtrim(result))
        output = result;
        return;
    end
end
end

function commands = get_nvidia_smi_commands()
if ~ispc
    commands = {'nvidia-smi'};
    return;
end

commands = {'nvidia-smi'};
candidatePaths = { ...
    fullfile(getenv('WINDIR'), 'System32', 'nvidia-smi.exe'), ...
    fullfile(getenv('ProgramFiles'), 'NVIDIA Corporation', 'NVSMI', 'nvidia-smi.exe'), ...
    fullfile(getenv('ProgramW6432'), 'NVIDIA Corporation', 'NVSMI', 'nvidia-smi.exe'), ...
    fullfile(getenv('ProgramFiles(x86)'), 'NVIDIA Corporation', 'NVSMI', 'nvidia-smi.exe') ...
    };

for iPath = 1:numel(candidatePaths)
    candidatePath = candidatePaths{iPath};
    if isempty(candidatePath) || exist(candidatePath, 'file') ~= 2
        continue;
    end
    quotedCommand = ['"' candidatePath '"'];
    if ~any(strcmp(commands, quotedCommand))
        commands{end + 1} = quotedCommand; %#ok<AGROW>
    end
end
end

function errorRedirect = get_shell_error_redirect()
if ispc
    errorRedirect = '2>nul';
else
    errorRedirect = '2>/dev/null';
end
end

function firstLine = get_first_line(textValue)
lineBreaks = strfind(textValue, newline);
if isempty(lineBreaks)
    firstLine = textValue;
else
    firstLine = textValue(1:lineBreaks(1)-1);
end
end

function cachedSelection = load_auto_device_cache()
cachedSelection = [];
cacheText = read_settings_cache();
if isempty(cacheText)
    cacheText = read_pref_cache();
end

if isempty(cacheText)
    return;
end

try
    cachedSelection = jsondecode(cacheText);
catch
    cachedSelection = [];
    return;
end

if ~is_valid_cache_for_this_pc(cachedSelection)
    cachedSelection = [];
end
end

function save_auto_device_cache(autoDeviceSelection)
if isempty(autoDeviceSelection) || ~isfield(autoDeviceSelection, 'usesAuto') || ~autoDeviceSelection.usesAuto
    return;
end

autoDeviceSelection.schemaVersion = 1;
autoDeviceSelection.machineId = get_machine_id();
autoDeviceSelection.platform = computer();
autoDeviceSelection.matlabRelease = version('-release');
autoDeviceSelection.savedAt = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss'));

cacheText = jsonencode(autoDeviceSelection);
if ~write_settings_cache(cacheText)
    setpref('xDDxAddon', 'AutoSimulationDevice', cacheText);
end
end

function clear_auto_device_cache()
clear_settings_cache();
if ispref('xDDxAddon', 'AutoSimulationDevice')
    rmpref('xDDxAddon', 'AutoSimulationDevice');
end
end

function tf = is_valid_cache_for_this_pc(cachedSelection)
tf = isstruct(cachedSelection) ...
    && isfield(cachedSelection, 'schemaVersion') ...
    && cachedSelection.schemaVersion == 1 ...
    && isfield(cachedSelection, 'machineId') ...
    && strcmp(cachedSelection.machineId, get_machine_id()) ...
    && isfield(cachedSelection, 'selectedDevice') ...
    && any(strcmp(cachedSelection.selectedDevice, {'cpu', 'cuda'}));
end

function cacheText = read_settings_cache()
try
    rootSettings = settings;
    cacheText = char(rootSettings.xDDxAddon.AutoSimulationDevice.ActiveValue);
catch
    cacheText = '';
end
end

function cacheWritten = write_settings_cache(cacheText)
try
    rootSettings = settings;
    try
        settingsGroup = rootSettings.xDDxAddon;
    catch
        addGroup(rootSettings, 'xDDxAddon');
        settingsGroup = rootSettings.xDDxAddon;
    end
    try
        settingsValue = settingsGroup.AutoSimulationDevice;
    catch
        addSetting(settingsGroup, 'AutoSimulationDevice', 'FactoryValue', '');
        settingsValue = settingsGroup.AutoSimulationDevice;
    end
    settingsValue.PersonalValue = cacheText;
    cacheWritten = true;
catch
    cacheWritten = false;
end
end

function clear_settings_cache()
try
    rootSettings = settings;
    rootSettings.xDDxAddon.AutoSimulationDevice.PersonalValue = '';
catch
end
end

function cacheText = read_pref_cache()
cacheText = '';
if ispref('xDDxAddon', 'AutoSimulationDevice')
    cacheText = getpref('xDDxAddon', 'AutoSimulationDevice');
end
end

function machineId = get_machine_id()
machineName = getenv('COMPUTERNAME');
if isempty(machineName)
    machineName = getenv('HOSTNAME');
end
if isempty(machineName)
    try
        machineName = char(java.net.InetAddress.getLocalHost.getHostName);
    catch
        machineName = 'unknown-host';
    end
end
machineId = [computer() ':' lower(strtrim(machineName))];
end
