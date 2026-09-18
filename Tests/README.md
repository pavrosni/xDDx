# xDDx tests

This suite mirrors the legacy `xDDx/Tests` intent for the reorganized tree and
adds deterministic unit and native CPU integration coverage. The normal suite
does not discover or execute the heterogeneous simulator. The explicitly
enabled external backend test runs its public `xDDx_simulator.m` entry script.

Run the normal suite from the project root:

```matlab
addpath('Tests')
results = run_all_tests
```

Collect statement coverage for the same normal suite and all eligible source
folders (requires MATLAB R2023a or newer):

```matlab
addpath('Tests')
[results, coverageTable] = run_test_coverage
```

The default excludes tests tagged `External`. To include Docker/backend and
legacy example comparisons:

```matlab
addpath('Tests')
results = run_all_tests('IncludeExternal', true)
```

For MATLAB R2017b, use the quoted name-value syntax shown above. The test
code avoids newer name-value call syntax, function-result dot indexing, and
multi-dimension reduction arguments. Compatibility changes are checked locally
on R2024a; execution on Linux/R2017b still needs verification, including the
library functions and external backends called by these tests.

Backend comparison environment variables:

- `XDDX_TEST_CPU_ARCHITECTURES`, default
  `sse2,avx,avx2,avx512,arm64`
- `XDDX_TEST_CUDA_VERSIONS`, default `cuda11,cuda12`
- `XDDX_TEST_EXECUTION_MODES`, default `native,docker`
- `XDDX_TEST_REPORT_DIRECTORY`, default
  `Tests/backend_comparison/reports`

The full backend test executes the exact numerical defaults in
`simulation_toolbox/heterogeneous_simulator/xDDx_simulator.m`. It only
disables the GUI and plots and explicitly selects the backend. Native cases
are Windows-only; Docker cases run wherever Docker is usable, including
Windows with Docker Desktop/WSL integration. Unsupported ISA, GPU, OS, and
Docker combinations are recorded as skipped with a reason.

Run the matrix directly:

```matlab
addpath('Tests/backend_comparison')
[backendResults, reportPaths, hostInfo] = run_xddx_backend_comparison
```

The runner updates timestamped CSV and MAT reports after every backend. Each
successful row contains the pressure-field L2 and L-infinity norms plus
absolute and relative L2/L-infinity differences from the first successful
backend. Reports also include timing, output size, host/OS/WSL/GPU/Docker
metadata, Git revision, failures, and skip reasons.

## Selecting simulation scripts

The detailed example-regression runner mirrors the legacy suite. Add its
folder to the path and list the available reorganized scripts:

```matlab
addpath('Tests')
list_simulation_scripts
```

Run one script against the legacy project:

```matlab
[summary, comparison, speed] = run_simulation_script_tests( ...
    'Scripts', 'fp_sf', ...
    'SimulationDevices', 'cpu');
```

By default, the reference project is the sibling `xDDx` checkout. Compare
the reorganized project directly with the oldest checkout used by the legacy
`xDDx/Tests` runner as follows:

```matlab
[summary, comparison, speed] = run_simulation_script_tests( ...
    'Scripts', 'fp_sf', ...
    'ReferenceProject', 'oldest', ...
    'SimulationDevices', 'cpu');
```

Reference choices are `xddx` (default), `oldest` (the sibling `xDDx_old`
checkout), and `custom`. For a custom checkout, provide both
`'ReferenceProject', 'custom'` and `'OldProjectDir', 'D:\path\to\checkout'`.
An explicitly supplied `OldProjectDir` remains backward-compatible and takes
precedence; known `xDDx` and `xDDx_old` paths are labeled automatically.

Run several scripts, with repeated timed runs:

```matlab
[summary, comparison, speed] = run_simulation_script_tests( ...
    'Scripts', {'quick_start_flat', 'transducer_simulation_sf'}, ...
    'SimulationDevices', {'cuda', 'cpu'}, ...
    'WarmupRuns', 1, ...
    'RunsPerScript', 3);
```

Run all mapped scripts using both CPU and CUDA variants:

```matlab
[summary, comparison, speed] = run_simulation_script_tests( ...
    'Scripts', 'all', ...
    'SimulationDevices', 'both');
```

`Scripts` accepts catalog names, `.m` filenames, reorganized relative paths,
or legacy relative paths. Use `'RunMode', 'new-only'` to execute only this
checkout. `NewExecutionMode` accepts `auto`, `native`, or `docker`.

The runner writes three timestamped CSV files matching the legacy report
families:

- `example_regression_summary_*.csv` - per script/device/project status,
  timing, Rayleigh core calls, outputs, and errors;
- `example_regression_comparison_*.csv` - sizes and numerical differences;
- `example_regression_speed_*.csv` - old/new timing ratios.

Every script under `holography_toolbox` and
`simulation_toolbox/homogeneous_simulatior` has an explicit legacy mapping.
For automated external test execution, `XDDX_EXAMPLE_REGRESSION_SCRIPTS`
still accepts a comma-separated selection. Set
`XDDX_EXAMPLE_REGRESSION_REFERENCE` to `xddx` or `oldest` to select its
reference project.

Reports are written only below `Tests/example_regression/reports`.
