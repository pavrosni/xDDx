% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% Returns an actionable hint when "docker run" failed, based on command output.
% Used by rayleigh_simulator to give the user guidance after a failed docker run.

function hint = get_docker_run_failure_hint(cmdout)

hint = '';
if isempty(cmdout)
    return;
end
outLower = lower(cmdout);
nl = char(10);

if isunix && ~ismac
    if ~isempty(strfind(outLower, 'permission denied')) || ...
       ~isempty(strfind(outLower, 'docker.sock')) || ...
       ~isempty(strfind(outLower, 'connect: permission denied'))
        hint = [ ...
            'Your user may need to be in the "docker" group to run without sudo:' nl ...
            '  1. Run in a terminal:  sudo usermod -aG docker $USER' nl ...
            '  2. Log out and log back in (or run:  newgrp docker)' nl ...
            '  3. Restart MATLAB and rerun the simulation.'];
        return;
    end
end

if ~isempty(strfind(outLower, 'cannot connect')) || ...
   ~isempty(strfind(outLower, 'daemon')) || ...
   ~isempty(strfind(outLower, 'is the docker daemon running'))
    hint = 'Start Docker: on Linux run "sudo systemctl start docker". On macOS start Docker Desktop. Then rerun the simulation.';
    return;
end

if ~isempty(strfind(outLower, 'no such image')) || ~isempty(strfind(outLower, 'pull')) || ~isempty(strfind(outLower, 'not found'))
    hint = 'Try pulling the image manually: docker pull <image>. Ensure your Docker username is set correctly in xddx_docker_config.m.';
end

end
