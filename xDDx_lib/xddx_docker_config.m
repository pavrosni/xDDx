% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% Configuration for xDDx Docker images used on Linux/Mac.
% Set dockerUsername to your Docker Hub username so the simulator can pull and run
% images: xddx-sf-cpu-sse2, xddx-sf-cpu-avx, xddx-sf-cpu-avx2,
% xddx-sf-cpu-avx512, xddx-sf-cpu-arm64, xddx-sf-cuda11, xddx-sf-cuda12,
% xddx-transient-cpu-sse2, xddx-transient-cpu-avx,
% xddx-transient-cpu-avx2, xddx-transient-cpu-avx512,
% xddx-transient-cpu-arm64, xddx-transient-cuda11,
% and xddx-transient-cuda12.

function [config] = xddx_docker_config

config = [];
config.dockerUsername = 'pavrosni';  
config.imageUpdatePeriodDays = 7;     

end
