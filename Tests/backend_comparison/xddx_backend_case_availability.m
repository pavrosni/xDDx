function [available, reason] = xddx_backend_case_availability( ...
    backendCase, hostInfo, projectRoot)
%XDDX_BACKEND_CASE_AVAILABILITY Check whether a case can run on this host.

if strcmp(backendCase.ExecutionMode, 'native') && ~ispc
    available = false;
    reason = 'Native simulator executables are Windows-only.';
    return;
end
if strcmp(backendCase.ExecutionMode, 'docker') ...
        && ~hostInfo.DockerAvailable
    available = false;
    reason = ['Docker is unavailable: ' hostInfo.DockerMessage];
    return;
end
if strcmp(backendCase.Device, 'cuda')
    [available, reason] = cuda_available(backendCase, hostInfo);
    if ~available
        return;
    end
else
    [available, reason] = cpu_variant_available(backendCase, hostInfo);
    if ~available
        return;
    end
end

if strcmp(backendCase.ExecutionMode, 'native')
    [available, reason] = native_files_available( ...
        backendCase, projectRoot);
else
    available = true;
end
end

function [available, reason] = cuda_available(backendCase, hostInfo)
available = false;
reason = '';
if ismac
    reason = 'CUDA is not supported on macOS.';
    return;
end
if ~hostInfo.HasCudaGpu
    reason = 'No NVIDIA GPU was detected with nvidia-smi.';
    return;
end

requiredVersion = str2double(backendCase.Variant(5:end));
if ~isnan(hostInfo.MaximumCudaVersion) ...
        && hostInfo.MaximumCudaVersion < requiredVersion
    reason = sprintf( ...
        'The NVIDIA driver reports CUDA %.1f, below %s.', ...
        hostInfo.MaximumCudaVersion, backendCase.Variant);
    return;
end
available = true;
end

function [available, reason] = cpu_variant_available(backendCase, hostInfo)
available = false;
reason = '';
if strcmp(backendCase.Variant, 'arm64')
    if ~hostInfo.IsArm64
        reason = 'ARM64 backend requested on a non-ARM64 host.';
        return;
    end
    available = true;
    return;
end
if hostInfo.IsArm64
    reason = 'x86 SIMD backend requested on an ARM64 host.';
    return;
end

if hostInfo.CpuFlagsKnown ...
        && ~cpu_flags_include_variant(hostInfo.CpuFlags, backendCase.Variant)
    reason = sprintf('CPU does not advertise %s support.', ...
        upper(backendCase.Variant));
    return;
end
available = true;
end

function tf = cpu_flags_include_variant(flags, variant)
switch variant
    case 'avx512'
        flagName = 'avx512f';
    otherwise
        flagName = variant;
end
tf = ~isempty(regexp([' ' lower(flags) ' '], ...
    ['\s' flagName '\s'], 'once'));
end

function [available, reason] = native_files_available( ...
    backendCase, projectRoot)
if strcmp(backendCase.Device, 'cpu')
    xddxExecutable = fullfile(projectRoot, 'xDDx_lib', ...
        'rayleigh_cpp', 'cpu', ...
        ['rayleigh_cpu-' backendCase.Variant '.exe']);
    kwaveExecutable = fullfile(projectRoot, 'simulation_toolbox', ...
        'heterogeneous_simulator', 'heterogeneous_core', 'k-Wave', ...
        'binaries', 'win_cpu', ...
        ['kspaceFirstOrder-OMP-' backendCase.Variant '.exe']);
else
    versionNumber = backendCase.Variant(5:end);
    xddxExecutable = fullfile(projectRoot, 'xDDx_lib', ...
        'rayleigh_cpp', 'cuda', ...
        ['rayleigh_cuda-' backendCase.Variant '.exe']);
    kwaveExecutable = fullfile(projectRoot, 'simulation_toolbox', ...
        'heterogeneous_simulator', 'heterogeneous_core', 'k-Wave', ...
        'binaries', 'win_cuda', ...
        ['kspaceFirstOrder-CUDA-' versionNumber '.exe']);
end

missingFiles = {};
if exist(xddxExecutable, 'file') ~= 2
    missingFiles{end + 1} = xddxExecutable;
end
if exist(kwaveExecutable, 'file') ~= 2
    missingFiles{end + 1} = kwaveExecutable;
end
available = isempty(missingFiles);
if available
    reason = '';
else
    reason = ['Missing native executable(s): ' ...
        strjoin(missingFiles, '; ')];
end
end
