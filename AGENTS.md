# Working with xDDx

This document provides instructions for AI agents working in the xDDx repository.

## Read first

- Read [README.md](README.md) for setup and the public entry point.
- Before configuring, running, or interpreting a heterogeneous simulation, read [docs/simulation-guide.md](docs/simulation-guide.md), then inspect the user's current `simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m` and relevant helpers.
- For Linux/macOS execution, also read [DOCKER_SETUP.md](DOCKER_SETUP.md). For code changes and test selection, read [Tests/README.md](Tests/README.md).
- Public repository URL: `https://github.com/pavrosni/xDDx`. Prefer the current local checkout and user-provided files over published defaults. If files or execution tools are unavailable, state what is missing; do not invent their contents or claim a run.

## Repository map

- `simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m`: user-editable entry script, inputs, execution, and plots.
- `simulation_toolbox/heterogeneous_simulator/lib/`: simulator core, input loaders, geometry checks, source generation, and visualization.
- `simulation_toolbox/heterogeneous_simulator/heterogeneous_core/k-Wave/`: integrated k-Wave MATLAB code and solver support.
- `xDDx_lib/`: xDDx propagation functions and backend support.
- `data_for_examples/`, `holography_toolbox/`, `simulation_toolbox/homogeneous_simulatior/`: data and related workflows; the last directory's spelling is intentional here.
- `Tests/`: unit, integration, and explicitly enabled external/backend tests.

## Working rules

- Users are expected to edit the simulator script. Inspect existing edits and preserve unrelated settings and work. For a simulation request, make focused input changes; change core algorithms only when the task requires it.
- Explain parameter changes using physical meaning, units, and effects. Ask for missing scientific inputs that materially affect the result; identify any assumptions.
- Run the script from `simulation_toolbox/heterogeneous_simulator` so its relative library paths resolve. Do not relocate a copy without adapting and checking its paths.
- Check MATLAB and the requested CPU/CUDA or Docker backend before execution. Start with one validated case before a large sweep; consider RAM, GPU memory, disk space, and runtime.
- Validate units, material dimensions, equal voxel spacing, target indices, boundary placement, and field limits. Preserve validation checks rather than bypassing failures.
- **Before relying on a new standard single-element simulation setup, perform `waterTest` and inspect the analytical comparison to check source geometry and field behavior.** Repeat after relevant geometry, source, frequency, or discretization changes. Follow the guide's apex and output-axis requirements. A run completing is not a validation pass.
- The built-in water test rejects arrays, custom sources, and untyped sources. Do not relabel them to bypass this restriction. Use an appropriate independent homogeneous-medium reference and geometry checks; report validation as incomplete when such a reference is unavailable.
- Save water-test evidence separately, restore the intended heterogeneous inputs, and rerun geometry checks before production. A water test does not validate the heterogeneous material map or establish convergence by itself.
- For unattended runs, explicitly control `useGUI` and `shouldPlot`. Disabling plots also disables the entry script's automatic water-test comparison; provide and inspect a separate comparison before claiming validation.
- `simulatorInputOverrides` permits only backend and UI fields listed in its helper. It cannot override `waterTest`, material, frequency, or geometry. Edit those inputs in the script; preassigning `simulatorInputs` before running it does not work because the script recreates that structure.
- Preserve input datasets and existing results. Save new runs to distinct locations with inputs, code revision/local edits, backend details, logs, and validation evidence. Do not mistake temporary solver files or workspace variables for an archived result.
- Report separately what was edited, executed, and validated, including warnings and limitations. Never report prepared inputs, skipped comparisons, cancelled runs, or empty results as successful simulations.

## Verification scope

Use checks appropriate to the change. Documentation-only changes need link/content checks, not solver runs. The normal test suite does not execute the heterogeneous simulator; passing it alone does not validate a simulator change. External backend comparisons are potentially expensive and are not a replacement for the water-test workflow.
