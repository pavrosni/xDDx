function [simulationDevice, autoDeviceSelection] = resolve_auto_simulation_device(simulationDevice, varargin)
%RESOLVE_AUTO_SIMULATION_DEVICE Resolve and cache the Rayleigh device flag.
%
% Use simulationDevice = 'auto' to select CUDA on CUDA-capable Windows/Linux
% machines and CPU otherwise. macOS is always resolved to CPU. The resolved
% choice, and the CUDA version when applicable, are cached in MATLAB settings
% after a successful simulation so later runs on the same PC avoid detection checks.

options = parse_options(varargin{:});

if options.SaveSelection
    save_auto_device_cache(options.AutoDeviceSelection);
    autoDeviceSelection = options.AutoDeviceSelection;
    return;
end

if options.ClearCache
    clear_auto_device_cache();
    autoDeviceSelection = [];
    return;
end

simulationDevice = normalize_simulation_device(simulationDevice);
autoDeviceSelection = create_empty_selection(simulationDevice);

if ~strcmp(simulationDevice, 'auto')
    return;
end

autoDeviceSelection.usesAuto = true;
cachedSelection = [];
if ~options.ForceDetect
    cachedSelection = load_auto_device_cache();
end

if ~isempty(cachedSelection)
    simulationDevice = cachedSelection.selectedDevice;
    autoDeviceSelection.usedCache = true;
    autoDeviceSelection.detectionSummary = cachedSelection.detectionSummary;
    autoDeviceSelection.selectedCudaVersion = get_cached_selection_cuda_version(cachedSelection);
else
    [simulationDevice, detectionSummary] = detect_auto_device();
    autoDeviceSelection.usedCache = false;
    autoDeviceSelection.detectionSummary = detectionSummary;
end

autoDeviceSelection.selectedDevice = simulationDevice;
if ~strcmp(simulationDevice, 'cuda')
    autoDeviceSelection.selectedCudaVersion = '';
end

if strcmp(simulationDevice, 'cuda') && ~isempty(autoDeviceSelection.selectedCudaVersion)
    fprintf('Auto simulation device: Rayleigh = %s, CUDA = %s (%s)\n', ...
        simulationDevice, autoDeviceSelection.selectedCudaVersion, autoDeviceSelection.detectionSummary);
else
    fprintf('Auto simulation device: Rayleigh = %s (%s)\n', ...
        simulationDevice, autoDeviceSelection.detectionSummary);
end
end

function options = parse_options(varargin)
options = struct();
options.ForceDetect = false;
options.SaveSelection = false;
options.ClearCache = false;
options.AutoDeviceSelection = [];

if mod(numel(varargin), 2) ~= 0
    error('resolve_auto_simulation_device:InvalidOptions', ...
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
            error('resolve_auto_simulation_device:InvalidOption', ...
                'Unknown option "%s".', optionName);
    end
end
end

function simulationDevice = normalize_simulation_device(simulationDevice)
if isa(simulationDevice, 'string')
    simulationDevice = char(simulationDevice);
end

if ~ischar(simulationDevice)
    error('resolve_auto_simulation_device:InvalidDevice', ...
        'simulationDevice must be ''auto'', ''cpu'', or ''cuda''.');
end

simulationDevice = lower(strtrim(simulationDevice));
if ~any(strcmp(simulationDevice, {'auto', 'cpu', 'cuda'}))
    error('resolve_auto_simulation_device:InvalidDevice', ...
        'simulationDevice must be ''auto'', ''cpu'', or ''cuda''.');
end
end

function selection = create_empty_selection(originalSimulationDevice)
selection = struct();
selection.schemaVersion = 2;
selection.usesAuto = false;
selection.usedCache = false;
selection.selectedDevice = originalSimulationDevice;
selection.selectedCudaVersion = '';
selection.originalSimulationDevice = originalSimulationDevice;
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

[status, cmdout] = system('nvidia-smi -L');
if status ~= 0
    return;
end

cmdoutLower = lower(cmdout);
if contains(cmdoutLower, 'gpu') || contains(cmdoutLower, 'nvidia')
    hasGpu = true;
    firstLine = get_first_line(strtrim(cmdout));
    summary = ['nvidia-smi detected CUDA-capable GPU: ' firstLine];
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
if isempty(autoDeviceSelection) ...
        || ~isfield(autoDeviceSelection, 'usesAuto') ...
        || ~autoDeviceSelection.usesAuto
    return;
end

autoDeviceSelection.schemaVersion = 2;
autoDeviceSelection.machineId = get_machine_id();
autoDeviceSelection.platform = computer();
autoDeviceSelection.matlabRelease = version('-release');
autoDeviceSelection.savedAt = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss'));
if ~isfield(autoDeviceSelection, 'selectedCudaVersion')
    autoDeviceSelection.selectedCudaVersion = '';
end

cacheText = jsonencode(autoDeviceSelection);
if ~write_settings_cache(cacheText)
    setpref('xDDxAddon', 'RayleighAutoSimulationDevice', cacheText);
end
end

function clear_auto_device_cache()
clear_settings_cache();
if ispref('xDDxAddon', 'RayleighAutoSimulationDevice')
    rmpref('xDDxAddon', 'RayleighAutoSimulationDevice');
end
end

function tf = is_valid_cache_for_this_pc(cachedSelection)
tf = isstruct(cachedSelection) ...
    && isfield(cachedSelection, 'schemaVersion') ...
    && any(cachedSelection.schemaVersion == [1 2]) ...
    && isfield(cachedSelection, 'machineId') ...
    && strcmp(cachedSelection.machineId, get_machine_id()) ...
    && isfield(cachedSelection, 'selectedDevice') ...
    && any(strcmp(cachedSelection.selectedDevice, {'cpu', 'cuda'})) ...
    && has_valid_cached_cuda_version(cachedSelection);
end

function cudaVersion = get_cached_selection_cuda_version(cachedSelection)
cudaVersion = '';
if isfield(cachedSelection, 'selectedCudaVersion') ...
        && is_valid_cuda_version(cachedSelection.selectedCudaVersion)
    cudaVersion = char(cachedSelection.selectedCudaVersion);
end
end

function tf = has_valid_cached_cuda_version(cachedSelection)
tf = true;
if ~isfield(cachedSelection, 'selectedCudaVersion') ...
        || isempty(cachedSelection.selectedCudaVersion)
    return;
end

tf = is_valid_cuda_version(cachedSelection.selectedCudaVersion);
end

function tf = is_valid_cuda_version(cudaVersion)
tf = any(strcmp(char(cudaVersion), {'cuda11', 'cuda12'}));
end

function cacheText = read_settings_cache()
try
    rootSettings = settings;
    cacheText = char(rootSettings.xDDxAddon.RayleighAutoSimulationDevice.ActiveValue);
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
        settingsValue = settingsGroup.RayleighAutoSimulationDevice;
    catch
        addSetting(settingsGroup, 'RayleighAutoSimulationDevice', 'FactoryValue', '');
        settingsValue = settingsGroup.RayleighAutoSimulationDevice;
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
    rootSettings.xDDxAddon.RayleighAutoSimulationDevice.PersonalValue = '';
catch
end
end

function cacheText = read_pref_cache()
cacheText = '';
if ispref('xDDxAddon', 'RayleighAutoSimulationDevice')
    cacheText = getpref('xDDxAddon', 'RayleighAutoSimulationDevice');
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
