function fix = ensure_xddx_docker_cli_in_path(candidateDockerPaths)
%ENSURE_XDDX_DOCKER_CLI_IN_PATH Add a discovered Docker CLI to MATLAB PATH.
%   macOS applications often inherit a shorter PATH than interactive
%   shells. This helper locates common Docker Desktop and Homebrew
%   installations and prepends the first matching directory.
%
%   candidateDockerPaths is optional and supports deterministic testing.

if nargin < 1
    candidateDockerPaths = { ...
        '/opt/homebrew/bin/docker', ...
        '/usr/local/bin/docker', ...
        '/Applications/Docker.app/Contents/Resources/bin/docker', ...
        '/usr/bin/docker'};
end

fix = struct('didChangePath', false, 'addedDirs', {{}});
currentPath = getenv('PATH');
pathDirectories = regexp(currentPath, pathsep, 'split');

for candidateIndex = 1:numel(candidateDockerPaths)
    dockerExecutable = char(candidateDockerPaths{candidateIndex});
    if exist(dockerExecutable, 'file') ~= 2
        continue;
    end

    dockerDirectory = fileparts(dockerExecutable);
    if any(strcmp(pathDirectories, dockerDirectory))
        continue;
    end

    if isempty(currentPath)
        updatedPath = dockerDirectory;
    else
        updatedPath = [dockerDirectory pathsep currentPath];
    end
    setenv('PATH', updatedPath);
    fix.didChangePath = true;
    fix.addedDirs = {dockerDirectory};
    return;
end
end
