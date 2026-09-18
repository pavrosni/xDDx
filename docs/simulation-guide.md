# xDDx simulation guide for users and AI assistants

Use this guide with the current [entry script](../simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m). Users may have changed its defaults. Read the actual script before proposing edits or executing it.

## Access and prerequisites

Public repository: `https://github.com/pavrosni/xDDx`. Keep the repository folder structure intact. A GitHub link supplies reference material, not access to a user's local edits or a MATLAB runtime. For chat assistance, provide this guide, [AGENTS.md](../AGENTS.md), the current script, and relevant helper files or error output. If an assistant cannot retrieve a file, it should request it explicitly.

MATLAB and working solver backends are needed for execution. Windows uses bundled native executables; Linux/macOS use [DOCKER_SETUP.md](../DOCKER_SETUP.md). Confirm the installed MATLAB release and required functions/products for the selected case; do not invent a minimum release or assume a GPU is available. An agent may use an available MATLAB integration or command-line execution. These Markdown files do not install or configure either.

## Choose and edit a case

1. Establish the medium, transducer, frequency, target, desired output region, and intended accuracy. Inspect existing local modifications.
2. Edit the relevant `MAIN PARAMETERS` and optional settings in the entry script. Keep unrelated settings intact. Inspect the input loaders and generator for accepted fields before constructing custom inputs.
3. Check geometry and resource requirements, then perform the water-test validation below before relying on a supported new setup.
4. Restore the heterogeneous configuration, run it, and save the inputs, results, and validation evidence together.

The entry script creates `simulatorInputs` afresh. Setting that variable in the workspace beforehand will not configure the run. Its optional `simulatorInputOverrides` accepts only `kWaveCalculationFlag`, `xDDxCalculationFlag`, `cpuArchitecture`, `cudaVersion`, `useGUI`, and `shouldPlot`; see [the override helper](../simulation_toolbox/heterogeneous_simulator/lib/apply_xddx_simulator_input_overrides.m). Physics inputs and `waterTest` must be changed in the script or supplied through a deliberately constructed complete call to the core. Do not assume the core supplies defaults omitted from the script.

## Inputs, units, and geometry

| Input | Meaning and checks |
| --- | --- |
| `MaterialMatrix` | Scalar struct with `c0`, `rho0`, `alpha`, `dx`, `dy`, `dz`. Property arrays have matching dimensions `[Nx, Ny, Nz]`, corresponding to x, y, z. |
| `c0`, `rho0` | Sound speed in m/s and density in kg/m^3. Check finite, physically meaningful values. |
| `alpha` | Attenuation in dB/cm at the transducer's operating frequency. The core performs the frequency-power conversion; do not apply it twice. Reassess attenuation when changing frequency. |
| `dx`, `dy`, `dz` | Positive voxel spacings in metres; the simulator requires `dx = dy = dz`. |
| `TransducerSf` | Source structure; frequency in Hz, dimensions in metres, velocity amplitudes in m/s, and array-element phases in radians. Use the generator/loader contract. |
| `ixTarget`, `iyTarget`, `izTarget` | One-based voxel indices, not distances. Check them against the material array dimensions. |
| `izBoundaryCondition` | Voxel index or `'auto'`. Check uniform contact-medium placement and its position relative to the source and target. Automatic selection differs between water and heterogeneous cases. |
| `xFieldBegin/End`, `yFieldBegin/End`, `zFieldBegin/End` | Output limits in metres, in the simulator's source coordinates. Equal begin/end requests a plane, line, or point, which must intersect actual grid nodes. |
| `CFL`, PML, source sampling | Numerical controls. Do not treat a successful run as evidence of spatial/time convergence or adequate absorbing boundaries. |

For transverse coordinates the core uses `x = (ix - ixTarget) * dx` and `y = (iy - iyTarget) * dy`. For a spherical source, `z = (iz - izTarget) * dz + radiusOfCurvature`; for a flat source, `z = (iz - izBoundaryCondition) * dz`. Field limits are therefore not simply distances from the first material voxel.

Consult [load_material_matrix.m](../simulation_toolbox/heterogeneous_simulator/lib/load_material_matrix.m), [generate_xDDx_transducer.m](../simulation_toolbox/heterogeneous_simulator/lib/generate_xDDx_transducer.m), and [load_xDDx_transducer.m](../simulation_toolbox/heterogeneous_simulator/lib/load_xDDx_transducer.m) for full input contracts. Keep runtime geometry and field-limit checks enabled.

## Required water-test validation

For a standard single-element source, run `simulatorInputs.waterTest = true` before relying on a new simulation setup. This replaces the medium inside the core with uniform contact-medium sound speed and density and zero attenuation. It helps check source geometry, coordinate placement, and axial field behavior against the analytical reference.

1. Preserve the intended heterogeneous configuration and note all temporary changes. Keep source dimensions, frequency, amplitudes, and discretization representative of the intended run.
2. Confirm `TransducerSf.type` is `standard_single_element`. The built-in test supports standard flat and spherical single-element sources. Arrays, custom sources, and untyped sources are rejected; do not change their type label to force acceptance.
3. For a spherical source, choose an in-bounds `izTarget` greater than `radiusOfCurvature / dz` so the bowl apex lies on the material grid. The entry script's current 50 mm radius and 0.5 mm spacing give a ratio of 100; its default `izTarget = 100` does not meet this documented condition. Do not simply toggle `waterTest` without checking placement. Record temporary shifts and recheck the original heterogeneous geometry when restoring it.
4. Set `waterTest = true` and normally use `izBoundaryCondition = 'auto'`. Include actual output grid points at `x = 0` and `y = 0`, and a nonzero axial interval (`zFieldBegin ~= zFieldEnd`). For a spherical source, retain an axial range reaching the focal region and satisfying the geometry checks.
5. Run with `shouldPlot = true` to obtain the entry script's automatic axial comparison. Use `useGUI = true` for an interactive setup review, or `false` for execution without the setup dialog. The comparison requires plotting even when the setup GUI is disabled.
6. Inspect the computed and analytical axial amplitude curves, focal position where applicable, and amplitude scale. Save the comparison and record discrepancies. The [comparison helper](../simulation_toolbox/heterogeneous_simulator/lib/plot_water_test_on_axis_comparison.m) draws a plot; it does not implement a numerical acceptance threshold. It warns and skips comparison if the output has no on-axis grid point. A missing plot, a skipped comparison, or solver completion alone is not a pass.
7. Investigate discrepancies in geometry, sampling, CFL, PML, and source parameters before relying on the result. State the acceptance criteria used; do not invent an established project tolerance. If quantitative or headless validation is needed, implement an explicit comparison using the helper's analytical reference and record its metrics and criteria. With `shouldPlot = false`, the script does not perform this comparison automatically.
8. Save water-test results separately. Restore `waterTest = false` and the intended material, target, boundary, output limits, and UI settings. Rerun geometry checks for the heterogeneous case. Water-test boundary selection and the resulting computational domain can differ from the production case.

For arrays and custom sources, explain that the built-in water test is unavailable. Check geometry and compare the actual source in a homogeneous medium against a suitable independent analytical/numerical reference. A single-element substitute does not validate the original array or custom geometry. If an appropriate reference is unavailable, report that validation remains incomplete.

Repeat the relevant validation after changing geometry, source, frequency, or discretization. Water validation does not establish correctness of the heterogeneous material map or full heterogeneous convergence.

## Execution

From the repository root in MATLAB:

```matlab
cd simulation_toolbox/heterogeneous_simulator
xDDx_simulator
```

Keep that working directory for the script's relative library paths. The default setup GUI waits for **Proceed**. In an unattended session, explicitly disable `useGUI`; disable `shouldPlot` only when the required validation comparison is handled separately. Check for stale workspace variables `simulatorInputOverrides` and `output_single_frequency_format`, which can change execution; clear or set these intentionally without clearing unrelated user work.

`'auto'` selects available backends. For reproducible comparisons, select explicit backends and record the actual CPU architecture/CUDA version or container image used. Estimate resources before increasing the grid or launching sweeps. Reducing the output dimensions does not necessarily reduce the computational domain proportionally.

`prepareSimulationOnly = true` requires `kWaveCalculationFlag = 'cpu'`. It prepares k-Wave input HDF5 and prints the external command, then returns without running that solver or producing a pressure field; earlier input/source preparation still consumes resources. Set `preparedSimulationDataPath` to a distinct output directory when preserving prepared files.

`showPresimulatedData` accepts an existing output HDF5 path in CPU mode. It cannot be combined with preparation mode. Preserve the matching inputs and geometry when loading saved output; an arbitrary output file is not interchangeable with a new configuration.

## Results and reporting

Check that `simulatorData` contains a nonempty, finite `complexPressure` field and sensor coordinates matching the requested region. Cancellation, some GUI validation failures, and preparation-only execution can return an empty structure. Inspect warnings and logs as well as numerical values.

`abs(simulatorData.complexPressure)` is pressure amplitude in Pa; `angle(...)` is phase in radians. These are workspace values; explicitly save results to archive them. Memory-saving mode can omit medium overlay data, so preserve original inputs independently.

Use a distinct run directory and record:

- The effective input configuration, source/material files or reproducible generation settings, and any random seeds.
- Repository revision and local changes, including a snapshot of the configured script.
- MATLAB release, selected backends, binary/container versions where available, and hardware.
- Pressure and sensor coordinates, logs, warnings, elapsed time, and output paths.
- Water-test or alternative validation evidence, criteria, discrepancies, and remaining limitations.

Summarize what changed, what actually executed, and what was validated. If MATLAB or a backend is unavailable, provide the configured files and reproducible commands and state that execution remains unverified.

## Code verification

See [Tests/README.md](../Tests/README.md) for commands and external/backend comparisons. The normal suite excludes heterogeneous simulator execution. Backend agreement is useful but does not replace geometry/water validation. For documentation-only changes, inspect links and consistency with source; no expensive solver run is needed.
