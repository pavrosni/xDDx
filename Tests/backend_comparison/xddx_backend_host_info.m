function hostInfo = xddx_backend_host_info(projectRoot)
%XDDX_BACKEND_HOST_INFO Collect machine metadata and backend capabilities.

hostInfo = struct();
hostInfo.MachineName = get_machine_name();
hostInfo.Platform = computer();
hostInfo.Architecture = computer('arch');
hostInfo.OperatingSystem = get_operating_system();
hostInfo.MatlabVersion = version();
hostInfo.MatlabRelease = version('-release');
hostInfo.CpuModel = get_cpu_model();
hostInfo.CpuFlags = get_cpu_flags();
hostInfo.CpuFlagsKnown = ~isempty(hostInfo.CpuFlags);
hostInfo.IsArm64 = detect_arm64();
hostInfo.IsWsl = detect_wsl();

[hostInfo.HasCudaGpu, hostInfo.NvidiaSummary, ...
    hostInfo.MaximumCudaVersion] = get_nvidia_info();
dockerCheck = check_docker_ready();
hostInfo.DockerAvailable = dockerCheck.ok;
hostInfo.DockerMessage = dockerCheck.message;
hostInfo.DockerHint = dockerCheck.hint;
hostInfo.DockerVersion = get_command_output('docker --version 2>&1');
hostInfo.GitCommit = get_git_commit(projectRoot);
hostInfo.CollectedUtc = utc_timestamp();
end

function value = get_machine_name()
value = getenv('COMPUTERNAME');
if isempty(value)
    value = getenv('HOSTNAME');
end
if isempty(value)
    try
        value = char(java.net.InetAddress.getLocalHost.getHostName);
    catch
        value = 'unknown-host';
    end
end
end

function value = get_operating_system()
try
    value = system_dependent('getos');
catch
    value = computer();
end
value = strtrim(value);
end

function value = get_cpu_model()
try
    value = char(feature('GetCPU'));
catch
    value = '';
end
if isempty(value)
    value = 'unknown';
end
end

function flags = get_cpu_flags()
flags = '';
if isunix && ~ismac
    try
        cpuText = fileread('/proc/cpuinfo');
        tokens = regexp(cpuText, ...
            '(?:flags|Features)\s*:\s*([^\n\r]+)', 'tokens', 'once');
        if ~isempty(tokens)
            flags = lower(strtrim(tokens{1}));
        end
    catch
        flags = '';
    end
elseif ismac
    flags = lower(get_command_output( ...
        'sysctl -n machdep.cpu.features machdep.cpu.leaf7_features 2>/dev/null'));
end
end

function tf = detect_arm64()
architectureText = lower([computer() ' ' computer('arch') ' ' ...
    getenv('PROCESSOR_ARCHITECTURE')]);
tf = contains(architectureText, 'arm64') ...
    || contains(architectureText, 'aarch64') ...
    || contains(architectureText, 'maca64');
if isunix
    machineArchitecture = lower(get_command_output('uname -m 2>/dev/null'));
    tf = tf || contains(machineArchitecture, 'arm64') ...
        || contains(machineArchitecture, 'aarch64');
end
end

function tf = detect_wsl()
tf = ~isempty(getenv('WSL_DISTRO_NAME')) ...
    || ~isempty(getenv('WSL_INTEROP'));
if ~tf && isunix && ~ismac
    kernelText = lower(get_command_output('uname -r 2>/dev/null'));
    tf = contains(kernelText, 'microsoft') || contains(kernelText, 'wsl');
end
end

function [hasGpu, summary, maximumCudaVersion] = get_nvidia_info()
hasGpu = false;
summary = get_command_output('nvidia-smi -L 2>&1');
maximumCudaVersion = NaN;
if isempty(summary) || contains(lower(summary), 'not recognized') ...
        || contains(lower(summary), 'not found')
    summary = '';
    return;
end

hasGpu = contains(lower(summary), 'gpu') ...
    || contains(lower(summary), 'nvidia');
fullOutput = get_command_output('nvidia-smi 2>&1');
tokens = regexp(fullOutput, ...
    'CUDA Version:\s*(\d+)\.(\d+)', 'tokens', 'once');
if ~isempty(tokens)
    maximumCudaVersion = str2double( ...
        [tokens{1} '.' tokens{2}]);
end
summary = first_line(summary);
end

function value = get_git_commit(projectRoot)
oldFolder = cd(projectRoot);
folderCleanup = onCleanup(@() cd(oldFolder));
value = get_command_output('git rev-parse HEAD 2>&1');
if contains(lower(value), 'fatal:')
    value = '';
end
end

function value = get_command_output(command)
[status, value] = system(command);
if status ~= 0
    value = '';
else
    value = strtrim(value);
end
end

function value = first_line(value)
lineBreak = regexp(value, '[\r\n]', 'once');
if ~isempty(lineBreak)
    value = value(1:lineBreak - 1);
end
end

function value = utc_timestamp()
timestamp = datetime('now', 'TimeZone', 'UTC');
timestamp.Format = 'yyyy-MM-dd''T''HH:mm:ssXXX';
value = char(timestamp);
end
