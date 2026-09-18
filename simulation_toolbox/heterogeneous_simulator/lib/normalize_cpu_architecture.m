function cpuArchitecture = normalize_cpu_architecture(cpuArchitecture)
%NORMALIZE_CPU_ARCHITECTURE Validate CPU architecture selection.

if nargin < 1 || isempty(cpuArchitecture)
    cpuArchitecture = 'auto';
end

if ~(ischar(cpuArchitecture) || (isstring(cpuArchitecture) && isscalar(cpuArchitecture)))
    error('cpuArchitecture must be ''auto'', ''avx512'', ''avx2'', ''avx'', ''sse2'', or ''arm64''.');
end

cpuArchitecture = lower(strtrim(char(cpuArchitecture)));
validArchitectures = {'auto', 'avx512', 'avx2', 'avx', 'sse2', 'arm64'};
if ~any(strcmp(cpuArchitecture, validArchitectures))
    error('cpuArchitecture must be ''auto'', ''avx512'', ''avx2'', ''avx'', ''sse2'', or ''arm64''.');
end

validate_cpu_architecture_platform(cpuArchitecture);

end

function validate_cpu_architecture_platform(cpuArchitecture)
if strcmp(cpuArchitecture, 'auto')
    return;
end

if strcmp(cpuArchitecture, 'arm64')
    if ismac && detect_apple_silicon()
        return;
    end
    if isunix && ~ismac && detect_linux_arm64()
        return;
    end
    if ispc && detect_windows_arm64()
        return;
    end
    error('cpuArchitecture = ''arm64'' is only supported on ARM64 platforms.');
end

if ismac && detect_apple_silicon()
    error('cpuArchitecture = ''%s'' is an x86/x64 variant and is not supported on Apple Silicon macOS. Use ''auto'' or ''arm64''.', cpuArchitecture);
end

if isunix && ~ismac && detect_linux_arm64()
    error('cpuArchitecture = ''%s'' is an x86/x64 variant and is not supported on ARM64 Linux. Use ''auto'' or ''arm64''.', cpuArchitecture);
end

end

function tf = detect_windows_arm64()
tf = ispc && (strcmpi(getenv('PROCESSOR_ARCHITECTURE'), 'ARM64') ...
    || strcmpi(getenv('PROCESSOR_ARCHITEW6432'), 'ARM64'));
end

function tf = detect_apple_silicon()
tf = false;

if ~ismac
    return;
end

[statusArch, archResult] = system('uname -m 2>/dev/null');
if statusArch == 0 && strcmp(strtrim(archResult), 'arm64')
    tf = true;
    return;
end

[statusBrand, brandResult] = system('sysctl -n machdep.cpu.brand_string 2>/dev/null');
if statusBrand == 0 && startsWith(strtrim(brandResult), 'Apple')
    tf = true;
end

end

function tf = detect_linux_arm64()
tf = false;

if ~(isunix && ~ismac)
    return;
end

[statusArch, archResult] = system('uname -m 2>/dev/null');
if statusArch == 0
    arch = lower(strtrim(archResult));
    tf = strcmp(arch, 'aarch64') || strcmp(arch, 'arm64');
end

end
