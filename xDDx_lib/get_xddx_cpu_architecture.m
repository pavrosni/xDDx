% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% Returns the CPU architecture suffix for xDDx binaries/images.
% ServiceParameters.cpuArchitecture may be 'auto', 'avx512', 'avx2',
% 'avx', 'sse2', or 'arm64'. Detection defaults to SSE2 if the CPU
% capabilities cannot be determined.

function [architecture, hasExplicitArchitecture] = get_xddx_cpu_architecture(ServiceParameters)

if nargin < 1
    ServiceParameters = [];
end

hasExplicitArchitecture = false;
if isstruct(ServiceParameters) && isfield(ServiceParameters, 'cpuArchitecture')
    requestedArchitecture = normalize_requested_cpu_architecture(ServiceParameters.cpuArchitecture);
    if ~strcmp(requestedArchitecture, 'auto')
        validate_requested_cpu_architecture_for_os(requestedArchitecture);
        architecture = requestedArchitecture;
        hasExplicitArchitecture = true;
        return;
    end
end

architecture = detect_cpu_architecture();
end

function architecture = normalize_requested_cpu_architecture(requestedValue)
if isa(requestedValue, 'string')
    requestedValue = char(requestedValue);
end

if ~ischar(requestedValue)
    error('ServiceParameters.cpuArchitecture must be ''auto'', ''avx512'', ''avx2'', ''avx'', ''sse2'', or ''arm64''.');
end

architecture = lower(strtrim(requestedValue));
validArchitectures = {'auto', 'avx512', 'avx2', 'avx', 'sse2', 'arm64'};
if ~any(strcmp(architecture, validArchitectures))
    error('ServiceParameters.cpuArchitecture must be ''auto'', ''avx512'', ''avx2'', ''avx'', ''sse2'', or ''arm64''.');
end
end

function validate_requested_cpu_architecture_for_os(architecture)
if ismac && ~strcmp(architecture, 'arm64')
    error('On macOS, ServiceParameters.cpuArchitecture must be ''auto'' or ''arm64''.');
end
end

function architecture = detect_cpu_architecture()
if ismac
    architecture = 'arm64';
    return;
elseif ispc
    architecture = 'avx2';
    return;
end

cpuFeatures = detect_cpu_features();
if cpuFeatures.avx512
    architecture = 'avx512';
elseif cpuFeatures.avx2
    architecture = 'avx2';
elseif cpuFeatures.avx
    architecture = 'avx';
else
    architecture = 'sse2';
end

end

function features = detect_cpu_features()
% Default to SSE2 if detection fails.
features = struct('sse2', false, 'avx', false, 'avx2', false, 'avx512', false);
if isunix && ~ismac
    % Linux: check /proc/cpuinfo flags
    features.sse2 = cpuinfo_has_flag('sse2');
    features.avx = cpuinfo_has_flag('avx');
    features.avx2 = cpuinfo_has_flag('avx2');
    features.avx512 = cpuinfo_has_flag('avx512f');
end
end

function hasFlag = cpuinfo_has_flag(flagName)
[status, ~] = system(['grep -q -w ' flagName ' /proc/cpuinfo 2>/dev/null']);
hasFlag = (status == 0);
end
