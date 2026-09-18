% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% Returns the Docker image name (without username) for Rayleigh simulation.
% Single-frequency images are named xddx-sf-cpu-sse2, xddx-sf-cpu-avx,
% xddx-sf-cpu-avx2, xddx-sf-cpu-avx512, xddx-sf-cpu-arm64,
% xddx-sf-cuda11, or xddx-sf-cuda12.
% Transient images are named xddx-transient-cpu-sse2,
% xddx-transient-cpu-avx, xddx-transient-cpu-avx2,
% xddx-transient-cpu-avx512, xddx-transient-cpu-arm64,
% xddx-transient-cuda11, or xddx-transient-cuda12.
% CPU detection defaults to SSE2, and CUDA detection defaults to CUDA12, if
% the respective capability cannot be determined.

function imageName = get_xddx_docker_image(simulationDevice, ServiceParameters, isTransient)

if nargin < 2
    ServiceParameters = [];
end
if nargin < 3
    isTransient = false;
end

if isTransient
    imagePrefix = 'xddx-transient';
else
    imagePrefix = 'xddx-sf';
end

if strcmp(simulationDevice, 'cpu')
    imageName = [imagePrefix '-cpu-' get_xddx_cpu_architecture(ServiceParameters)];
elseif strcmp(simulationDevice, 'cuda')
    imageName = [imagePrefix '-' get_xddx_cuda_version(ServiceParameters)];
else
    error('simulationDevice must be ''cpu'' or ''cuda''.');
end

end
