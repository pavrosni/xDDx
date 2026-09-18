# Backend selection and CUDA troubleshooting

Use this guide for slow automatic detection, CUDA failures, and NVIDIA driver updates. Settings below belong in the existing input block of [xDDx_simulator.m](../simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m).

## Slow first-run backend detection

With backend settings left at `'auto'`, the first run may spend several minutes detecting CPU/CUDA support and selecting an executable version. This delay has been observed with older MATLAB releases and NVIDIA drivers; MATLAB may show **Busy** during detection. The backend choice is cached after a successful simulation, so subsequent runs on the same computer normally avoid the initial device-detection delay.

If MATLAB appears stuck during detection, first try **Ctrl+C**. If it remains unresponsive, terminate the MATLAB process (for example, using Windows Task Manager) and restart MATLAB. Unsaved workspace data and editor changes may be lost. Before running again, edit the existing backend settings in [xDDx_simulator.m](../simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m) to select the device and CPU architecture or CUDA version explicitly.

If you are unsure which settings to choose, use this conservative CPU configuration on **Windows/Linux computers with Intel or AMD x86-64 processors**:

```matlab
simulatorInputs.kWaveCalculationFlag = 'cpu';
simulatorInputs.xDDxCalculationFlag = 'cpu';
simulatorInputs.cpuArchitecture = 'sse2';
```

These settings bypass automatic CPU/CUDA and CPU architecture selection for both cores. `cudaVersion` is unused in CPU mode. Change these assignments in the script itself, including `cpuArchitecture` in the optional settings; assigning them in the workspace before running the script will not override its settings.

For explicit GPU selection, set both calculation flags to `'cuda'` and set `simulatorInputs.cudaVersion` to `'cuda11'` or `'cuda12'`, matching an available executable and your NVIDIA driver's capabilities. For the macOS CPU configuration, follow the [Docker setup instructions](DOCKER_SETUP.md).

## Check and update NVIDIA drivers on Windows

If CUDA fails, check the NVIDIA graphics driver first. MATLAB detecting a GPU does not guarantee that the driver supports the supplied CUDA executables.

### 1. Detect the driver

Open Windows Command Prompt or PowerShell and run:

```text
nvidia-smi
```

Look for **Driver Version** and **CUDA Version**. The CUDA version here describes what the driver supports; it is not the installed CUDA Toolkit version. See [NVIDIA's documentation](https://docs.nvidia.com/deploy/nvidia-smi/index.html).

### 2. Check CUDA support

| Reported CUDA Version | What to do |
| --- | --- |
| 12.x or newer | Try `'cuda12'` first. If it fails, try `'cuda11'`. |
| 11.x | Use `'cuda11'`, or update the driver to use CUDA 12. |
| Below 11 | Update the driver before using either supplied CUDA version, or use CPU mode. |
| Command not found, no GPU listed, or driver communication error | The driver may be missing or not working. Check **Device Manager > Display adapters** for an NVIDIA GPU, then install or update its driver. A missing command alone does not prove the driver is absent. |

The GPU must also be supported by the selected executable. Even when `nvidia-smi` reports CUDA 12 or newer, the CUDA 11 build may work where the CUDA 12 build fails because the builds can differ in GPU architecture support and runtime requirements. Keep both calculation flags set to `'cuda'`, change `simulatorInputs.cudaVersion` to `'cuda11'`, and rerun. If both versions fail, use the CPU settings above and keep the full MATLAB error message for troubleshooting.

### 3. Install or update the driver

Go to the [official NVIDIA driver download page](https://www.nvidia.com/en-us/drivers/), select your GPU model and Windows version, and download the driver offered for that combination. For laptops that require a manufacturer-specific driver, use the laptop manufacturer's official support page.

Save your work, close MATLAB, run the installer, and restart Windows. Install the **graphics driver**.

### 4. Run again

Run `nvidia-smi` again to confirm the driver update and its CUDA support. In [xDDx_simulator.m](../simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m), set:

```matlab
simulatorInputs.kWaveCalculationFlag = 'cuda';
simulatorInputs.xDDxCalculationFlag = 'cuda';
simulatorInputs.cudaVersion = 'cuda12'; % Use 'cuda11' for CUDA 11.x support.
```

Run the script again from `simulation_toolbox/heterogeneous_simulator`. Check that it finishes without errors and produces a pressure field; reaching 100% alone does not establish success.

For Linux and Docker setup, see [DOCKER_SETUP.md](DOCKER_SETUP.md). macOS uses the CPU backend.
