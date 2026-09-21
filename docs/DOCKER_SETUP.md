# Linux and macOS setup

**Linux and macOS:** install Docker first. Docker runs the precompiled solvers in packaged environments called *containers*, while you continue working in MATLAB as usual.

## macOS

Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/) for your chip (Apple silicon or Intel), open it from Applications, and leave it running. Optionally, in Docker Desktop, open [Settings > General](https://docs.docker.com/desktop/settings-and-maintenance/settings/#general) and enable **Start Docker Desktop when you sign in to your computer** so it starts automatically at login.

## Linux

Follow the [Docker Engine installation guide](https://docs.docker.com/engine/install/) for your distribution, start the Docker service, and complete the [non-root setup](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user) so your user can run Docker without `sudo`. Log out and back in after changing group membership.

## Check and run on either platform

1. **Restart MATLAB and check access.** On macOS, ensure that Docker Desktop is running.
2. **Run an example of your choice**, such as [xDDx_simulator.m](../simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m) or [quick_start_spherical.m](../holography_toolbox/quick_start_spherical.m). Before running it, navigate to the corresponding folder in MATLAB. Stay connected to the internet during the first run: the toolbox automatically downloads the required xDDx and k-Wave solver images and reuses them on later runs. Keep the default Docker image settings.

If a Docker container fails to start, see the official troubleshooting guides for [Docker Desktop on macOS](https://docs.docker.com/desktop/troubleshoot-and-support/troubleshoot/) or [Docker Engine on Linux](https://docs.docker.com/engine/daemon/troubleshoot/).

**Optional Linux GPU support:** install an NVIDIA driver and [configure the NVIDIA Container Toolkit for Docker](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html). Until GPU access is configured, set both `simulatorInputs.xDDxCalculationFlag` and `simulatorInputs.kWaveCalculationFlag` to `'cpu'`. macOS uses the CPU backend.

## Control Docker updates

By default, each solver checks for image updates when a simulation runs and **7 days** have passed since its last update attempt. If an update fails, an existing local image is used.

To change this, edit `config.imageUpdatePeriodDays` in **both** [xddx_docker_config.m](xDDx_lib/xddx_docker_config.m) and [kwave_docker_config.m](simulation_toolbox/heterogeneous_simulator/heterogeneous_core/k-Wave/kwave_docker_config.m):

```matlab
config.imageUpdatePeriodDays = Inf; % Disable recurring update checks
```

Use a positive number such as `30` for a 30-day interval, or restore `7` for the default. `Inf` keeps using downloaded images after the initial update attempt; missing images or missing update history still trigger a download attempt.

**To update manually**, open your system terminal and run:

```sh
docker image ls
docker pull REPOSITORY:TAG
```

Replace `REPOSITORY:TAG` with the repository and tag shown for each solver image you use, and repeat for both xDDx and k-Wave. You can find the full list of repositories on the [Docker Hub page for xDDx and Memory-Saving k-Wave](https://hub.docker.com/u/pavrosni). Choose the appropriate variants for your hardware: CPU (SSE2, AVX, AVX2, AVX-512), ARM64, or CUDA (CUDA 11, CUDA 12). For example, for the AVX2 CPU versions:
```sh
docker pull pavrosni/xddx-sf-cpu-avx2:latest
docker pull pavrosni/kspacefirstorder-omp-avx2:latest
```
The next simulation will use the updated images, even if recurring checks are disabled. See the [Docker pull reference](https://docs.docker.com/reference/cli/docker/image/pull/) for details.
