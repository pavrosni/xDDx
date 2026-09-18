%%PATHS TO LIBRARIES AND BINARIES

%    Library folders relative to the MATLAB current folder. Run this script
%    from simulation_toolbox/heterogeneous_simulator so these paths resolve.
%    xDDx library for acoustic field propagation and source boundary conditions.
xDDxLibPath = '../../xDDx_lib';
%    Integrated k-Wave library for wave propagation through heterogeneous media.
kWavePath = 'heterogeneous_core/k-Wave';
%    Heterogeneous simulator helpers for input loading, geometry checks, 
%    and visualization.
xDDxSimulatorLibPath = 'lib';

simulatorInputs = struct();
simulatorInputs.scriptDirectory = fileparts(mfilename('fullpath'));

%    Add paths to the MATLAB path
addpath(genpath(kWavePath));
addpath(genpath(xDDxLibPath));
addpath(genpath(xDDxSimulatorLibPath));

simulatorInputs.kWaveBinPath = initialize_simulator_paths(simulatorInputs.scriptDirectory, kWavePath);

%% INPUT BLOCK

% =========================================================================
% MAIN PARAMETERS
% Start here: these settings define the basic simulation setup most users
% need to choose first. The optional sections below are still important, but
% their defaults are suitable for a first run.
% =========================================================================

% 1) Calculation flags
%    Use 'auto' to select CUDA when available and CPU otherwise.
%    Use 'cuda' or 'cpu' to force a specific backend.
simulatorInputs.kWaveCalculationFlag = 'auto';
simulatorInputs.xDDxCalculationFlag = 'auto';

% 2) Medium parameters
%    The medium is defined on a Cartesian grid: z is the transducer axis,
%    and x and y are the transverse directions.
% 
%   'MaterialMatrix' scalar struct with fields:
%       'c0': 3D matrix with the sound speed in m/s at each grid node
%       of the heterogeneous medium
%       'rho0': 3D matrix with the density in kg/m^3 at each grid node
%       of the heterogeneous medium
%       'alpha': 3D matrix with the attenuation in dB/cm at the operating
%       frequency of the transducer (TransducerSf.frequency, given in Hz)
%       at each grid node of the heterogeneous medium. The simulator
%       converts these values to power law coefficients by dividing by
%       (TransducerSf.frequency * 1e-6)^alphaPower, where alphaPower is
%       specified separately in the simulator inputs
%       'dx': x-step of the medium Cartesian grid in m (positive scalar)
%       'dy': y-step of the medium Cartesian grid in m (positive scalar)
%       'dz': z-step of the medium Cartesian grid in m (positive scalar)
%
%    All fields are required and must be non-empty and numeric. The matrices
%    'c0', 'rho0', and 'alpha' must have identical sizes [Nx, Ny, Nz], with
%    the first, second, and third dimensions corresponding to x, y, and z,
%    respectively (ndgrid format).
% 
%    The current source model requires dx = dy = dz. Reinterpolate data
%    onto a grid with equal spacing before using it here.

%    Option A (default): build material matrix in memory
%    (an automatically generated skull model)
    simulatorInputs.MaterialMatrix = generate_skull_material_matrix('dx', 0.5e-3,'dy', 0.5e-3,'dz', 0.5e-3);

%    Option B: load from MAT
%    (file must contain variable MaterialMatrix, see load_material_matrix)
%     simulatorInputs.MaterialMatrix = load_material_matrix(materialMatrixPath);

%    Contact medium: the homogeneous coupling medium between the transducer
%    and the heterogeneous body (typically water).
%    Match them to the coupling region in MaterialMatrix. 
%    Sound speed in m/s, 1500 m/s is a typical water approximation.
%    Density in kg/m^3, 1000 kg/m^3 is a typical water approximation.
simulatorInputs.soundSpeedContactMedium = 1500;
simulatorInputs.densityContactMedium = 1000;

% 3) CFL number
%    Courant-Friedrichs-Lewy number: CFL = c_max * dt / dx, the ratio of
%    the distance a wave travels in one time step to the grid spacing.
%    It is a dimensionless time step. Here, the maximum sound speed in
%    the medium is used. 0.3 is a good starting point.
%    Reduce it to test numerical stability.
%    With the new memory-saving k-Wave cores, lowering CFL does not
%    increase memory consumption.
simulatorInputs.CFL = 0.3;

% 4) Target and boundary condition parameters
%    x index of the target
simulatorInputs.ixTarget = 217;
%    y index of the target
simulatorInputs.iyTarget = 217;
%    z index of the target. For a spherical water-test, this must be
%    greater than radiusOfCurvature/dz so the bowl apex lies on-grid.
simulatorInputs.izTarget = 100;
%    z index of the boundary condition. Use 'auto' to select the last
%    uniform slice before the first non-uniform slice along z.
simulatorInputs.izBoundaryCondition = 'auto';

% 5) Transducer parameters
%    Choose one source for TransducerSf:
%
%    Option A (default): build a single-element spherical transducer.
    simulatorInputs.TransducerSf = generate_xDDx_transducer(...
        'aperture', 50e-3, ...
        'radiusOfCurvature', 50e-3, ...
        'frequency', 1e6, ...
        'initialVelocity', 1/simulatorInputs.soundSpeedContactMedium/simulatorInputs.densityContactMedium);

%    Option B: build a multi-element spherical array with circular elements.
%    This example uses the same 50 mm aperture, 50 mm radius of curvature,
%    and 1 MHz frequency as Option A. Twenty 6 mm circular elements are placed
%    at randomized positions. All elements have equal velocity and zero phase. 
%   simulatorInputs.TransducerSf = generate_xDDx_transducer(...
%       'aperture', 50e-3, ...
%       'radiusOfCurvature', 50e-3, ...
%       'frequency', 1e6, ...
%       'elementAperture', 6e-3, ...
%       'xCenters', [1, 7, -15, -14, -13, -5, 13, 1, 16, -6, 14, 9, -4, -5, 20, -10, 5, 9, -19, 2] * 1e-3, ...
%       'yCenters', [-1, -15, 1, -10, 11, -10, 3, 13, -11, 19, 13, -4, -18, 4, 0, -5, 6, 18, -4, -8] * 1e-3, ...
%       'initialVelocityAmplitudes', ones(20, 1)/simulatorInputs.soundSpeedContactMedium/simulatorInputs.densityContactMedium, ...
%       'initialVelocityPhases', zeros(20, 1));

%    Option C: load custom transducer from MAT
%    (file must contain variable TransducerSf; see load_xDDx_transducer)
%    transducerMatPath = 'my_transducer_sf.mat';
%     simulatorInputs.TransducerSf = load_xDDx_transducer(transducerMatPath);

% 6) 3D output window
%    x, y, and z limits in m for the rectangular output region.
%    Set "Begin" equal to "End" for a specific coordinate to reduce
%    the number of output dimensions.
simulatorInputs.xFieldBegin = -30e-3;
simulatorInputs.xFieldEnd   =  30e-3;
simulatorInputs.yFieldBegin = -30e-3; 
simulatorInputs.yFieldEnd   =  30e-3;
simulatorInputs.zFieldBegin =  10e-3;
simulatorInputs.zFieldEnd   =  70e-3;

% =========================================================================
% OPTIONAL PARAMETERS
% The settings below tune plotting, memory behavior, PML, technical 
% visualization details, CPU architecture, and CUDA version. 
% Adjust them for specific studies, or leave the defaults unchanged while 
% getting started.
% =========================================================================

% Visualization parameters
%    Use GUI for error handling and simulation setup validation
simulatorInputs.useGUI = true; 
%    If true, plot figures in this script after core computation
simulatorInputs.shouldPlot = true;

% Operation regimes
%    If true, drop medium overlay data after setup to minimize RAM.
simulatorInputs.strongMemorySavingMode = true;
%    CPU mode only: write input HDF5, print the external run command,
%    and stop before launching the solver.
simulatorInputs.prepareSimulationOnly = false;
%    CPU mode only: full path to an existing k-Wave *_output.h5 file to
%    load and plot instead of running the solver.
simulatorInputs.showPresimulatedData = false;
%    Folder for prepared input/output files when prepareSimulationOnly is
%    true. Leave empty to use the default k-Wave location.
simulatorInputs.preparedSimulationDataPath = '';
%    If true, substitute MaterialMatrix with a uniform medium:
%    MaterialMatrix.c0 = soundSpeedContactMedium,
%    MaterialMatrix.rho0 = densityContactMedium,
%    MaterialMatrix.alpha = 0.
%    Then compare k-Wave with the on-axis analytical solution (O'Neil).
simulatorInputs.waterTest = false;


%PML parameters
simulatorInputs.xSizePML  = 10;  % in number of grid points
simulatorInputs.ySizePML  = 10;  % in number of grid points
simulatorInputs.zSizePML  = 10;  % in number of grid points
simulatorInputs.xAlphaPML = 2;   % in Nepers/grid point
simulatorInputs.yAlphaPML = 2;   % in Nepers/grid point
simulatorInputs.zAlphaPML = 2;   % in Nepers/grid point
simulatorInputs.alphaPower = 2;  % power of the alpha coefficient 


%Technical parameters
%    Isolevels related to the pressure maximum for extracting
%    the isosurface of the simulated 3D field.
levelArray = [0.8 0.5 0.3];
transparencyArray = [0.6 0.5 0.2]; % transparency of the isosurfaces for the levels from levelArray
pressureTransparency = 0.5; % transparency of the pressure overlay
%    Radial reserves for the boundary condition width in x and y:
%    Wx -> (1 + radialReserveX) * Wx.
%    Wy -> (1 + radialReserveY) * Wy.
simulatorInputs.radialReserveX = 0.2;
simulatorInputs.radialReserveY = 0.2;
%    Minimum back-projection distance in wavelengths for the flat
%    boundary condition.
simulatorInputs.shiftHoloDistInWl = 10;
simulatorInputs.cpuArchitecture = 'auto'; % 'auto', 'avx512', 'avx2', 'avx', 'sse2', or 'arm64'
                                          % Used by CPU xDDx and k-Wave binaries.
simulatorInputs.cudaVersion = 'auto'; % 'auto', 'cuda11', 'cuda12', 11, or 12
                                      % Used by CUDA xDDx and k-Wave binaries.


%% END OF INPUT BLOCK

% A caller may provide a narrow set of overrides while still executing this
% exact default case. This is primarily used by the backend comparison test
% to disable interactive output and select one explicit backend.
if exist('simulatorInputOverrides', 'var')
    simulatorInputs = apply_xddx_simulator_input_overrides( ...
        simulatorInputs, simulatorInputOverrides);
end


%% MAIN BLOCK
% Technical checks and preparations before running the simulator core.
if ~isequal(simulatorInputs.showPresimulatedData, false)
    simulatorInputs.useGUI = false;
end

if exist('output_single_frequency_format', 'var')
    simulatorInputs.output_single_frequency_format = output_single_frequency_format;
end

% Run the simulator core
simulatorData = xDDx_simulator_core(simulatorInputs);

% Visualize the results when plotting is enabled and core returned plot data.
if simulatorInputs.shouldPlot && ~isempty(fieldnames(simulatorData))
    % Visualize the vibrational velocity at the boundary
    figure;
    imagesc(simulatorData.yGridSimulationVec([1 end])*1e3, simulatorData.xGridSimulationVec([1 end])*1e3, abs(simulatorData.VelocityHologram));
    xlabel('y, mm');
    ylabel('x, mm');
    axis equal;
    axis tight;
    colormap jet;
    colorbar;
    title(['Vibrational velocity in m/s at the boundary z = ' num2str_significant_figures(simulatorData.zGridSimulationVec(1) * 1e3, simulatorData.dx * 1e3) ' mm']);
    
    % Visualize the field
    visualize_field(simulatorData.xSensor, simulatorData.ySensor, simulatorData.zSensor, simulatorData.complexPressure, simulatorData.MaterialMatrixForPlot, ...
        simulatorData.xGridMediumVecForPlot, simulatorData.yGridMediumVecForPlot, simulatorData.zGridMediumVecForPlot, simulatorInputs.ixTarget, simulatorInputs.iyTarget, ...
        simulatorData.xBoundaryCondition, simulatorData.yBoundaryCondition, simulatorData.zBoundaryCondition, ...
        'levelArray', levelArray, ...
        'transparencyArray', transparencyArray, ...
        'pressureTransparency', pressureTransparency, ...
        'aperture', [simulatorData.apertureNominalX simulatorData.apertureNominalY], ...
        'radiusOfCurvature', simulatorData.radiusOfCurvature, ...
        'xDDxPath', xDDxLibPath, ...
        'transducerSf', simulatorData.TransducerSf);
    
    % For water test, plot the on-axis comparison with the analytical solution.
    if simulatorInputs.waterTest && simulatorInputs.zFieldBegin ~= simulatorInputs.zFieldEnd
        plot_water_test_on_axis_comparison(simulatorData.complexPressure, ...
           simulatorData.xGridSimulationVec, simulatorData.yGridSimulationVec, simulatorData.zGridSimulationVec, ...
           simulatorInputs.xFieldBegin, simulatorInputs.xFieldEnd, simulatorInputs.yFieldBegin, simulatorInputs.yFieldEnd, simulatorInputs.zFieldBegin, simulatorInputs.zFieldEnd, ...
           simulatorData.TransducerSf, simulatorInputs.soundSpeedContactMedium, simulatorInputs.densityContactMedium, ...
           simulatorData.frequency, simulatorData.aperture, simulatorData.isSphericalSource, simulatorData.radiusOfCurvature, simulatorData.dz, simulatorInputs.izTarget);
    end
end
    
