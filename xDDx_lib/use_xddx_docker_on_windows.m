% Copyright (c) 2026, the code is written by Pavel Rosnitskiy

function on = use_xddx_docker_on_windows()
%USE_XDDX_DOCKER_ON_WINDOWS Return whether Windows Docker mode is enabled.
%   XDDX_USE_DOCKER_ON_WINDOWS is the canonical environment variable.
%   XDDX_USE_DOCKER remains supported as a backward-compatible fallback
%   when the canonical variable is unset. The canonical variable takes
%   precedence when both are set.

environmentValue = strtrim(getenv('XDDX_USE_DOCKER_ON_WINDOWS'));
if isempty(environmentValue)
    environmentValue = strtrim(getenv('XDDX_USE_DOCKER'));
end

enabledValues = {'1', 'true', 'yes'};
on = ispc && ismember(lower(environmentValue), enabledValues);

end
