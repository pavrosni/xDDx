# xDDx: Ultrasound Simulation and Acoustic Holography

For a quick start, all toolbox files are stored in xDDx.zip, available on the releases page in the Assets drop-down list below the release notes.

xDDx is a numerical toolbox for ultrasound field simulation, transducer characterization, and acoustic holography. Its main entry point is the **heterogeneous simulator**, which combines xDDx source projection with a new, memory-saving integrated [k-Wave solver](#k-wave-license). The toolbox also provides acoustic holography and homogeneous-medium simulation tools.

- [xDDx Heterogeneous Simulator](#xddx-heterogeneous-simulator)
- [xDDx Acoustic Holography and Homogeneous Simulator](#xddx-acoustic-holography-and-homogeneous-simulator)
- [Documentation and Help](#documentation-and-help)
- [Citation and License](#citation-and-license)

## xDDx Heterogeneous Simulator

**Simulate ultrasound fields through heterogeneous media, including voxelized body models, from a single MATLAB script.** [xDDx_simulator.m](simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m) computes single-frequency pressure amplitude and phase for your transducer and a voxelized map of sound speed, density, and attenuation.

### Why use it

- **Reduce memory use.** For spherically focused transducers, an xDDx-powered Adaptive Boundary Condition (ABC) provides a flat boundary condition on the native k-Wave Cartesian grid, reducing the simulation volume. An integrated memory-saving k-Wave core optimized for single-frequency calculations further reduces memory use.
- **Model your source.** Generate flat or spherical single-element transducers and arrays with individual element amplitudes and phases, or load a custom `TransducerSf`, including a surface velocity distribution reconstructed by acoustic holography.
- **Run on CPU or GPU.** Both xDDx and k-Wave use compiled C++ backends with CPU and CUDA GPU versions. Automatic selection uses CUDA when available and CPU otherwise.
- **Inspect the results.** Preview the simulation geometry, then explore pressure isosurfaces and interactive slices. Save amplitude and phase for your own analysis.

### How it works

For a spherical transducer, xDDx projects its field onto a flat ABC plane in the contact medium, just before the heterogeneous region. This represents the curved source without including the whole transducer-to-body gap in the k-Wave grid. The integrated k-Wave solver propagates the wave through your material model and returns the complex pressure in your chosen output region. Its new CPU/GPU core generates harmonic source signals from compact amplitude and phase data and supports direct amplitude/phase output, reducing memory use and output file size for single-frequency simulations.

### Setup

**Windows:** with MATLAB installed, the complete toolbox is ready to use with its bundled executables. GPU execution also requires an NVIDIA driver compatible with the supplied CUDA executables, see [driver installation and verification](docs/BACKEND_SELECTION.md#check-and-update-nvidia-drivers-on-windows).

**Linux and macOS:** see the [Docker setup and update instructions](docs/DOCKER_SETUP.md).

#### Slow first run?
On the first run, automatic detection may leave MATLAB showing **Busy** for several minutes. This is an initial setup delay: after a successful simulation, the backend choice is cached, so subsequent runs on the same computer normally skip detection. See [first-run troubleshooting](docs/BACKEND_SELECTION.md#slow-first-run-backend-detection) for help with detection delays and explicit backend settings.

#### CUDA failure?
The CUDA backend may fail if your drivers are not up to date. A detected GPU or a progress counter reaching 100% does not guarantee a successful CUDA calculation. See [backend selection and troubleshooting](docs/BACKEND_SELECTION.md) for driver checks, installation steps, and explicit CPU/CUDA settings.

### Quick Start

Download **xDDx.zip** from the **Assets** section of the releases page. Unpack the archive into a folder of your choice, keeping its directory structure intact.

In MATLAB, navigate to the `xDDx/simulation_toolbox/heterogeneous_simulator` folder, open [xDDx_simulator.m](simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m), and click **Run** in the MATLAB Editor. Once setup is complete, the included example is ready to run out of the box.

Comments in `xDDx_simulator.m` guide you through the script and explain its settings if you want to explore the details or customize the simulation.

The defaults generate a synthetic skull model and a **1 MHz spherical transducer with a 50 mm aperture and 50 mm radius of curvature**. No input dataset is needed. Review the setup window and click **Proceed** to calculate and plot the field.

The results are returned in the MATLAB workspace as `simulatorData` and visualized automatically. Save them explicitly if you want to keep them after the session:

```matlab
pressureAmplitude = abs(simulatorData.complexPressure); % Pa
pressurePhase = angle(simulatorData.complexPressure);   % radians
```

### Make it your own

Edit the **MAIN PARAMETERS** in [xDDx_simulator.m](simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m), then rerun:

- **Medium:** replace the generated skull with your own `MaterialMatrix`. Use equally spaced voxels (`dx = dy = dz`), supply attenuation in **dB/cm at the operating frequency**. See the [input format](simulation_toolbox/heterogeneous_simulator/lib/load_material_matrix.m).
- **Transducer:** try the commented array example or load a custom source. See [transducer options](xDDx_lib/generate_xDDx_transducer.m).
- **Target and output:** choose target voxel indices and field limits in metres. Request a volume, plane, line, or point.

Before relying on a new standard single-element setup, follow the [water-test validation workflow](docs/simulation-guide.md#required-water-test-validation). Setting `simulatorInputs.waterTest = true` replaces the body model with a uniform medium for comparison with O'Neil's analytical solution. Check the source placement and output axis before running, and inspect the analytical comparison; enabling the flag alone is not validation. Arrays and custom sources require an independent homogeneous-medium reference.

The optional settings also support memory-saving operation, CPU input preparation for external execution, and reloading saved solver output.

### xDDx Simulator is AI-friendly

You are welcome to use an AI chat or coding agent to help you explore the simulator, choose parameters, or run simulations. The xDDx toolbox includes built-in [agent instructions](AGENTS.md) and a [simulation guide](docs/simulation-guide.md) to help your AI assistant follow the project workflow. Open this repository in your coding agent, or share the repository link and your current script with an AI chat. Try this prompt:

> Help me configure xDDx for [describe your simulation]. Read AGENTS.md and docs/simulation-guide.md in this repository (or at `https://github.com/pavrosni/xDDx`) first, and use my current xDDx_simulator.m as the starting point. Explain the parameter changes, validate the geometry with waterTest where supported, and tell me what you actually ran and checked. If you cannot access a required file, ask me for it.

## xDDx Acoustic Holography and Homogeneous Simulator

### Overview

The acoustic holography and homogeneous simulation tools provide transducer characterization and field projection through a scripting interface and interactive graphical tools. These tools were developed for MATLAB or Octave on Windows computers.

- **Acoustic holography:** forward or backward projection of transient acoustic holography measurements to identify transducer surface defects and reveal structural details of the radiated acoustic field. An automated procedure corrects axis misalignments to improve visualization of transducer surface vibrations.
- **Homogeneous simulation:** simulation of fields radiated by user-defined transducer designs in a uniform medium without attenuation, using the same projection algorithm as the holography tools.

The core algorithm is based on the Rayleigh integral implemented in C++ executables, with versions for CUDA-compatible GPUs and CPUs. The algorithms are optimized for post-processing planar scan data into holograms and for field projection calculations.

### Installation and Tool Organization

Use the toolbox archive and setup instructions in the [heterogeneous simulator section](#setup). In this repository, the tools are organized as follows:

| Tools | Folder |
| --- | --- |
| Acoustic holography | [holography_toolbox](holography_toolbox/) |
| Homogeneous simulation | [simulation_toolbox/homogeneous_simulatior](simulation_toolbox/homogeneous_simulatior/) |

Holography tools are further grouped by task, for example in [simple_projection_tools](holography_toolbox/simple_projection_tools/).

Open the corresponding `.m` script and run it from its folder in MATLAB/Octave. GUI tools assist with some interactive steps. The script comments describe each tool's capabilities, inputs, and outputs, and the included examples provide starting configurations. When adapting an example, keep an original copy and preserve the relative paths to libraries and data.

### Quick Start

Start with these two examples for an overview of the interface, inputs, and outputs:

1. **Acoustic holography:** [quick_start_spherical.m](holography_toolbox/quick_start_spherical.m) demonstrates hologram alignment and back-projection.
2. **Homogeneous simulation:** [transducer_simulation_sf.m](simulation_toolbox/homogeneous_simulatior/transducer_simulation_sf.m) calculates the field of a user-defined transducer at a single frequency.

## Documentation and Help

- **Heterogeneous simulation:** read the comments in [xDDx_simulator.m](simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m) and consult the [xDDx User Manual](https://github.com/pavrosni/xDDx/releases/download/v.1.0.0-alpha/xDDx_user_manual.pdf). For additional guidance on inputs, geometry checks, validation, and saving results, see the optional [simulation guide](docs/simulation-guide.md).
- **Acoustic holography and homogeneous simulation:** see the [xDDx User Manual](https://github.com/pavrosni/xDDx/releases/download/v.1.0.0-alpha/xDDx_user_manual.pdf), also available in the release Assets and toolbox archive, for detailed workflows and capabilities.
- **Linux and macOS setup:** see [Docker setup and update instructions](docs/DOCKER_SETUP.md).
- **Development and verification:** see [Tests/README.md](Tests/README.md).

For questions about the manual or toolbox, or to report bugs, contact Pavel Rosnitskiy at [pavrosni@gmail.com](mailto:pavrosni@gmail.com).

## Citation and License

If you find the **xDDx Heterogeneous Simulator** useful for your work, please consider citing the following paper:

> P. B. Rosnitskiy, O. A. Sapozhnikov, T. D. Khokhlova, and V. A. Khokhlova, "Memory-saving version of k-Wave toolbox for single-frequency simulations," _IEEE Trans. Ultrasonics_, 2026 (Early Access).

If you find **xDDx Acoustic Holography and Homogeneous Simulator** useful for your work, please consider citing the following paper:

> P. B. Rosnitskiy, O. A. Sapozhnikov, V. A. Khokhlova, W. Kreider, S. A. Tsysar, G. P. L. Thomas, K. Contreras, and T. D. Khokhlova, "xDDx: a Numerical Toolbox for Ultrasound Transducer Characterization and Design with Acoustic Holography," _IEEE Trans. Ultrason., Ferroelectr., Freq. Control_, vol. 72, no. 5, pp. 564–580, May 2025. DOI: [10.1109/TUFFC.2025.3542405](https://doi.org/10.1109/TUFFC.2025.3542405).

and the acoustic holography paper:

> O. A. Sapozhnikov, S. A. Tsysar, V. A. Khokhlova, and W. Kreider, "Acoustic holography as a metrological tool for characterizing medical ultrasound sources and fields," _The Journal of the Acoustical Society of America_, vol. 138, no. 3. Acoustical Society of America (ASA), pp. 1515–1532, Sep. 01, 2015.

Consult the toolbox distribution's `license.txt` for general license information.

### k-Wave License

The heterogeneous simulator includes a modified k-Wave solver maintained in the [pavrosni/k-wave fork](https://github.com/pavrosni/k-wave). k-Wave is distributed under the GNU Lesser General Public License (LGPL). See the bundled [k-Wave license information and recommended citations](simulation_toolbox/heterogeneous_simulator/heterogeneous_core/k-Wave/helpfiles/k-wave_license.html) for details.
