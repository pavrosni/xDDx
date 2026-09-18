% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% Returns the CUDA executable/image suffix supported by the current machine:
% 'cuda11' or 'cuda12'. Detection uses the maximum CUDA runtime version
% reported by nvidia-smi, which reflects the installed NVIDIA driver.
% If detection is unavailable, default to CUDA 12 to match the Docker path
% except on Windows 7/older, where CUDA 12 is not a practical target.

function [cudaVersion, hasExplicitCudaVersion] = get_xddx_cuda_version(ServiceParameters)

if nargin < 1
    ServiceParameters = [];
end

if ismac
    error('CUDA simulation is not supported on macOS. Use CPU mode instead (simulationDevice = ''cpu'').');
end

hasExplicitCudaVersion = false;
if isstruct(ServiceParameters) && isfield(ServiceParameters, 'cudaVersion')
    requestedCudaVersion = normalize_requested_cuda_version(ServiceParameters.cudaVersion);
    if ~strcmp(requestedCudaVersion, 'auto')
        cudaVersion = requestedCudaVersion;
        hasExplicitCudaVersion = true;
        return;
    end
end

cudaVersion = detect_cuda_version();
end

function cudaVersion = normalize_requested_cuda_version(requestedValue)
if isnumeric(requestedValue)
    requestedValue = num2str(requestedValue);
elseif isa(requestedValue, 'string')
    requestedValue = char(requestedValue);
end

if ~ischar(requestedValue)
    error('ServiceParameters.cudaVersion must be ''auto'', ''cuda11'', ''cuda12'', 11, or 12.');
end

requestedValue = lower(strtrim(requestedValue));
if strcmp(requestedValue, '11')
    requestedValue = 'cuda11';
elseif strcmp(requestedValue, '12')
    requestedValue = 'cuda12';
end

validCudaVersions = {'auto', 'cuda11', 'cuda12'};
if ~any(strcmp(requestedValue, validCudaVersions))
    error('ServiceParameters.cudaVersion must be ''auto'', ''cuda11'', ''cuda12'', 11, or 12.');
end

cudaVersion = requestedValue;
end

function cudaVersion = detect_cuda_version()
cudaVersion = get_default_cuda_version();

nvidiaSmiOutput = get_nvidia_smi_output();
if isempty(nvidiaSmiOutput)
    return;
end

% Parse "CUDA Version: X.Y" from nvidia-smi output.
cudaTokens = regexp(nvidiaSmiOutput, 'CUDA Version:\s*(\d+)\.(\d+)', 'tokens', 'once');
if ~isempty(cudaTokens)
    majorVersion = str2double(cudaTokens{1});
    if ~isnan(majorVersion) && majorVersion < 12
        cudaVersion = 'cuda11';
    elseif ~isnan(majorVersion)
        cudaVersion = 'cuda12';
    end
    return;
end

% Older nvidia-smi versions may not print "CUDA Version". Fall back to the
% NVIDIA driver version: CUDA 12 on Windows requires R527 or newer.
driverTokens = regexp(nvidiaSmiOutput, 'Driver Version:\s*(\d+)', 'tokens', 'once');
if isempty(driverTokens)
    return;
end

driverMajorVersion = str2double(driverTokens{1});
if ~isnan(driverMajorVersion) && driverMajorVersion < 527
    cudaVersion = 'cuda11';
elseif ~isnan(driverMajorVersion)
    cudaVersion = 'cuda12';
end
end

function cudaVersion = get_default_cuda_version()
if ispc && is_windows_7_or_older()
    cudaVersion = 'cuda11';
else
    cudaVersion = 'cuda12';
end
end

function output = get_nvidia_smi_output()
output = '';
commands = get_nvidia_smi_commands();

for iCommand = 1:numel(commands)
    [status, result] = system([commands{iCommand} ' ' get_shell_error_redirect()]);
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

function tf = is_windows_7_or_older()
tf = false;
if ~ispc
    return;
end

[status, result] = system('ver');
if status ~= 0
    return;
end

versionTokens = regexp(result, 'Version\s+(\d+)\.(\d+)', 'tokens', 'once');
if isempty(versionTokens)
    versionTokens = regexp(result, '\[(?:Version\s+)?(\d+)\.(\d+)', 'tokens', 'once');
end
if isempty(versionTokens)
    return;
end

majorVersion = str2double(versionTokens{1});
minorVersion = str2double(versionTokens{2});
tf = ~isnan(majorVersion) && ~isnan(minorVersion) ...
    && (majorVersion < 6 || (majorVersion == 6 && minorVersion <= 1));
end

function errorRedirect = get_shell_error_redirect()
if ispc
    errorRedirect = '2>nul';
else
    errorRedirect = '2>/dev/null';
end
end
