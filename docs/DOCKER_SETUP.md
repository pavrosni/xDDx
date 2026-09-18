# Linux and macOS setup

**Linux and macOS:** install Docker first. Docker runs the precompiled solvers in packaged environments called *containers*, while you continue working in MATLAB as usual.

## macOS

Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/) for your chip (Apple silicon or Intel), open it from Applications, and leave it running. Optionally, in Docker Desktop, open [Settings > General](https://docs.docker.com/desktop/settings-and-maintenance/settings/#general) and enable **Start Docker Desktop when you sign in to your computer** so it starts automatically at login.

## Linux

Follow the [Docker Engine installation guide](https://docs.docker.com/engine/install/) for your distribution, start the Docker service, and complete the [non-root setup](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user) so your user can run Docker without `sudo`. Log out and back in after changing group membership.

## Check and run on either platform

1. **Restart MATLAB and check access.** Run this in the MATLAB Command Window:

   ```matlab
   system('docker ps')
   ```

   A container list (even an empty one) and a return value of `0` mean Docker is ready.
2. **Run the [included example](README.md#try-the-included-example).** Keep an internet connection for the first run: the toolbox automatically downloads the required xDDx and k-Wave solver images and reuses them on later runs. Keep the default Docker image settings.

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

Replace `REPOSITORY:TAG` with the repository and tag shown for each solver image you use, and repeat for both xDDx and k-Wave. You can find the full list of repositories on [Docker Hub](https://hub.docker.com/u/pavrosni). For example, for the AVX2 CPU versions:

```sh
docker pull pavrosni/xddx-sf-cpu-avx2:latest
docker pull pavrosni/kspacefirstorder-omp-avx2:latest
```

Use your own listed variants for ARM64 or CUDA. The next simulation uses the updated images, even with recurring checks disabled. See the [Docker pull reference](https://docs.docker.com/reference/cli/docker/image/pull/) for details.
