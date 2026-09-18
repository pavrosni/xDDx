% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% Pre-flight check for Docker on Linux/Mac. Returns struct with:
%   ok       - true if Docker is installed and usable (without sudo on Linux)
%   message  - short error message if ok is false
%   hint     - actionable hint for the user (how to fix the problem)

function [result] = check_docker_ready

result = [];
result.ok = true;
result.message = '';
result.hint = '';

% 1) Check if Docker is installed
[status, ~] = system('docker --version 2>&1');
if status ~= 0 && (ismac || isunix)
    pathFix = ensure_xddx_docker_cli_in_path();
    if pathFix.didChangePath
        [status, ~] = system('docker --version 2>&1');
    end
end
if status ~= 0
    result.ok = false;
    result.message = 'Docker is not installed or not found in PATH.';
    result.hint = [ ...
        'Install Docker and ensure the "docker" command is available in the system PATH. ' ...
        'Linux: install Docker Engine (https://docs.docker.com/engine/install/). ' ...
        'macOS: install Docker Desktop (https://docs.docker.com/desktop/install/mac-install/). ' ...
        'Then restart MATLAB and rerun the simulation.'];
    return;
end

% 2) Check if we can run Docker (daemon running, and on Linux: no sudo required)
[status, cmdout] = system('docker ps 2>&1');
if status ~= 0
    result.ok = false;
    outLower = lower(cmdout);

    % Linux: permission denied / docker.sock -> suggest adding user to docker group
    if isunix && ~ismac
        if ~isempty(strfind(outLower, 'permission denied')) || ...
           ~isempty(strfind(outLower, 'docker.sock')) || ...
           ~isempty(strfind(outLower, 'connect: permission denied'))
            result.message = 'MATLAB cannot run Docker without sudo.';
            nl = char(10);
            result.hint = [ ...
                'Add your user to the "docker" group so you can run Docker without sudo:' nl ...
                '  1. Run in a terminal:  sudo usermod -aG docker $USER' nl ...
                '  2. Log out and log back in (or run:  newgrp docker)' nl ...
                '  3. Restart MATLAB, then rerun the simulation.' nl ...
                'If you run MATLAB from a launcher, start it after logging in again so it sees the new group.'];
            return;
        end
    end

    % Docker daemon not running (Linux or Mac)
    if ~isempty(strfind(outLower, 'cannot connect')) || ...
       ~isempty(strfind(outLower, 'daemon')) || ...
       ~isempty(strfind(outLower, 'is the docker daemon running'))
        result.message = 'Docker daemon is not running.';
        result.hint = [ ...
            'Start Docker: on Linux run "sudo systemctl start docker" (or start Docker Desktop if you use it). ' ...
            'On macOS start Docker Desktop from Applications. Then rerun the simulation.'];
        return;
    end

    % Generic failure
    result.message = 'Docker is installed but "docker ps" failed.';
    result.hint = 'Check that Docker Desktop (Mac) or the Docker service (Linux) is running. Then rerun the simulation.';
end

end
