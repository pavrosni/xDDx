function simulatorData = xDDx_simulator_core(simulatorInputs)
% Run the non-visual simulator workflow and return data for plotting.

% The current source model requires equal grid spacing in all directions.
dx = simulatorInputs.MaterialMatrix.dx;
dy = simulatorInputs.MaterialMatrix.dy;
dz = simulatorInputs.MaterialMatrix.dz;
if dx ~= dy || dx ~= dz
    error('xDDx:UnequalGridSpacing', ...
        ['The current source model requires equal grid steps (dx = dy = dz). ', ...
        'Received dx = %.16g, dy = %.16g, dz = %.16g m. ', ...
        'Reinterpolate the material data onto a grid with equal spacing in all three directions before running the simulation.'], ...
        dx, dy, dz);
end

originalSimulatorInputs = simulatorInputs;
[simulatorInputs, autoDeviceSelection] = resolve_auto_simulation_devices(simulatorInputs);

try
    simulatorData = xDDx_simulator_core_run(simulatorInputs);
    if shouldSaveAutoSelection(autoDeviceSelection, simulatorData)
        resolve_auto_simulation_devices(simulatorInputs, ...
            'AutoDeviceSelection', autoDeviceSelection, ...
            'SaveSelection', true);
    end
catch ME
    if autoDeviceSelection.usesAuto && autoDeviceSelection.usedCache
        resolve_auto_simulation_devices(simulatorInputs, ...
            'AutoDeviceSelection', autoDeviceSelection, ...
            'ClearCache', true);
        [simulatorInputs, retryAutoDeviceSelection] = resolve_auto_simulation_devices(originalSimulatorInputs, ...
            'ForceDetect', true);
        try
            simulatorData = xDDx_simulator_core_run(simulatorInputs);
            if shouldSaveAutoSelection(retryAutoDeviceSelection, simulatorData)
                resolve_auto_simulation_devices(simulatorInputs, ...
                    'AutoDeviceSelection', retryAutoDeviceSelection, ...
                    'SaveSelection', true);
            end
            return;
        catch retryME
            throwAutoSelectionError(retryME);
        end
    elseif autoDeviceSelection.usesAuto
        throwAutoSelectionError(ME);
    else
        rethrow(ME);
    end
end
end

function simulatorData = xDDx_simulator_core_run(simulatorInputs)
% Run the simulator after manual or automatic device flags are resolved.

kWaveBinPath = simulatorInputs.kWaveBinPath;
scriptDirectory = simulatorInputs.scriptDirectory;
kWaveCalculationFlag = simulatorInputs.kWaveCalculationFlag;
xDDxCalculationFlag = simulatorInputs.xDDxCalculationFlag;
MaterialMatrix = simulatorInputs.MaterialMatrix;
soundSpeedContactMedium = simulatorInputs.soundSpeedContactMedium;
densityContactMedium = simulatorInputs.densityContactMedium;
ixTarget = simulatorInputs.ixTarget;
iyTarget = simulatorInputs.iyTarget;
izTarget = simulatorInputs.izTarget;
izBoundaryCondition = simulatorInputs.izBoundaryCondition;
TransducerSf = simulatorInputs.TransducerSf;
validate_TransducerSf_type(TransducerSf);
xSizePML = simulatorInputs.xSizePML;
ySizePML = simulatorInputs.ySizePML;
zSizePML = simulatorInputs.zSizePML;
xAlphaPML = simulatorInputs.xAlphaPML;
yAlphaPML = simulatorInputs.yAlphaPML;
zAlphaPML = simulatorInputs.zAlphaPML;
CFL = simulatorInputs.CFL;
alphaPower = simulatorInputs.alphaPower;
xFieldBegin = simulatorInputs.xFieldBegin;
xFieldEnd = simulatorInputs.xFieldEnd;
yFieldBegin = simulatorInputs.yFieldBegin;
yFieldEnd = simulatorInputs.yFieldEnd;
zFieldBegin = simulatorInputs.zFieldBegin;
zFieldEnd = simulatorInputs.zFieldEnd;
radialReserveX = simulatorInputs.radialReserveX;
radialReserveY = simulatorInputs.radialReserveY;
shiftHoloDistInWl = simulatorInputs.shiftHoloDistInWl;
useGUI = simulatorInputs.useGUI;
strongMemorySavingMode = simulatorInputs.strongMemorySavingMode;
prepareSimulationOnly = simulatorInputs.prepareSimulationOnly;
showPresimulatedData = simulatorInputs.showPresimulatedData;
preparedSimulationDataPath = simulatorInputs.preparedSimulationDataPath;
waterTest = simulatorInputs.waterTest;

if isfield(simulatorInputs, 'cpuArchitecture')
    cpuArchitecture = normalize_cpu_architecture(simulatorInputs.cpuArchitecture);
else
    cpuArchitecture = 'auto';
end

if isfield(simulatorInputs, 'cudaVersion')
    cudaVersion = normalize_cuda_version(simulatorInputs.cudaVersion);
else
    cudaVersion = 'auto';
end

ServiceParameters = [];
ServiceParameters.threadsPerBlockGPU = 128;
if ~strcmp(cpuArchitecture, 'auto')
    ServiceParameters.cpuArchitecture = cpuArchitecture;
end
if ~strcmp(cudaVersion, 'auto')
    ServiceParameters.cudaVersion = cudaVersion;
end

outputSingleFrequencyFormat = [];
if isfield(simulatorInputs, 'output_single_frequency_format')
    outputSingleFrequencyFormat = simulatorInputs.output_single_frequency_format;
end

simulatorData = struct();

% Extract transducer parameters
if isfield(TransducerSf, 'radiusOfCurvature')
    radiusOfCurvature = TransducerSf.radiusOfCurvature;
else
    radiusOfCurvature = [];
end
frequency = TransducerSf.frequency;
[apertureNominalX, apertureNominalY] = estimate_source_aperture_nominal(TransducerSf);
aperture = max(apertureNominalX, apertureNominalY);
apertureXY = [apertureNominalX apertureNominalY];
isSphericalSource = ~isempty(radiusOfCurvature);

if waterTest
    validate_water_test_transducer(TransducerSf);
    MaterialMatrix.c0 = soundSpeedContactMedium * ones(size(MaterialMatrix.c0));
    MaterialMatrix.rho0 = densityContactMedium * ones(size(MaterialMatrix.rho0));
    MaterialMatrix.alpha = zeros(size(MaterialMatrix.alpha));
end

% Set up the binary path for k-Wave based on the k-Wave calculation flag
if ispc
    kWaveBinPath = fullfile(kWaveBinPath, ['win_', lower(kWaveCalculationFlag)]);
else
    kWaveBinPath = fullfile(kWaveBinPath, ['docker_', lower(kWaveCalculationFlag)]);
end

[prepareSimulationOnly, presimulatedOutputFile, preparedSimulationDataPath] = ...
    plan_presimulation_workflow(kWaveCalculationFlag, prepareSimulationOnly, showPresimulatedData, ...
    preparedSimulationDataPath, scriptDirectory);

% Resolve k-Wave binary/archive name from the CPU architecture and CUDA version flags.
if isfield(simulatorInputs, 'kWaveBinName') && ~strcmpi(simulatorInputs.kWaveBinName, 'auto')
    kWaveBinName = simulatorInputs.kWaveBinName;
else
    kWaveBinName = getKWaveBinaryName(kWaveCalculationFlag, kWaveBinPath, cpuArchitecture, cudaVersion);
end

% Auto-select boundary condition position (flat, spherical + heterogeneous, or spherical + water-test apex)
if ischar(izBoundaryCondition) && strcmpi(izBoundaryCondition, 'auto')
    izBoundaryCondition = select_boundary_condition_auto(MaterialMatrix, ...
        'IsSpherical', isSphericalSource, ...
        'WaterTest', waterTest, ...
        'IzTarget', izTarget, ...
        'RadiusOfCurvature', radiusOfCurvature);
    fprintf('izBoundaryCondition (auto): %d\n', izBoundaryCondition);
end

if ~validate_boundary_condition_geometry( ...
        MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, ...
        apertureXY, radiusOfCurvature, radialReserveX, radialReserveY, ...
        xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, useGUI)
    return;
end

% Validate the simulation window limits.
validate_field_limits(xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd);

% Define the medium grid: number of grid points in each direction and grid step.
nxMedium = size(MaterialMatrix.c0,1);
nyMedium = size(MaterialMatrix.c0,2);
nzMedium = size(MaterialMatrix.c0,3);
dx = MaterialMatrix.dx;
dy = MaterialMatrix.dy;
dz = MaterialMatrix.dz;

% Build shifted medium grid vectors directly from indices.
xGridMediumVec = ((1:nxMedium).' - ixTarget) * dx;
yGridMediumVec = ((1:nyMedium).' - iyTarget) * dy;
if isSphericalSource
    zGridMediumVec = ((1:nzMedium).' - izTarget) * dz + radiusOfCurvature;
else
    zGridMediumVec = ((1:nzMedium).' - izBoundaryCondition) * dz;
end

if isSphericalSource
    izApex = izTarget - (radiusOfCurvature/dz);
    if mod(izApex, 2) ~= 0
        izApex = fix(izApex) + 1;
    end
end

% Generate the simulation model: pad or crop the material matrix to the size of the simulation grid.
alphaContactMedium = 0; % in dB/cm (xDDx doesn't work with absorptive medium)

MaterialMatrixSimulation = generate_simulation_model(radialReserveX, radialReserveY, apertureXY, radiusOfCurvature, ...
    zGridMediumVec, izBoundaryCondition, dx, dy, dz, izTarget, zFieldEnd, ...
    xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, ...
    xSizePML, ySizePML, zSizePML, ixTarget, iyTarget, ...
    nxMedium, nyMedium, nzMedium, MaterialMatrix, ...
    soundSpeedContactMedium, densityContactMedium, alphaContactMedium, strongMemorySavingMode);

% Convert dB/cm to dB/cm/MHz.
MaterialMatrixSimulation.alpha = MaterialMatrixSimulation.alpha/(frequency*1e-6)^alphaPower;

% Define the simulation grid.
nxSimulationGrid = size(MaterialMatrixSimulation.c0,1);
nySimulationGrid = size(MaterialMatrixSimulation.c0,2);
nzSimulationGrid = size(MaterialMatrixSimulation.c0,3);

zBoundaryOffset = zGridMediumVec(izBoundaryCondition);
[keepOverlayData, clearMediumDataEarly] = plan_visualization_overlay_data(strongMemorySavingMode, useGUI);

if clearMediumDataEarly
    clear MaterialMatrix xGridMediumVec yGridMediumVec zGridMediumVec
end

% Shift the simulation grid to have the boundary condition coordinates at the z-origin.
kGridSimulation = kWaveGrid(nxSimulationGrid, dx, nySimulationGrid, dy, nzSimulationGrid, dz);
xGridSimulationVec = kGridSimulation.x_vec;
yGridSimulationVec = kGridSimulation.y_vec;
zGridSimulationVec = kGridSimulation.z_vec - kGridSimulation.z_vec(1) + zBoundaryOffset;

[xBoundaryCondition, yBoundaryCondition] = ndgrid(xGridSimulationVec, yGridSimulationVec);
zBoundaryCondition = zGridSimulationVec(1) + zeros(nxSimulationGrid, nySimulationGrid);

if useGUI
    % Let user validate target/boundary/box geometry before running simulation.
    proceed = validate_simulation_setup_window(MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, ...
        apertureXY, radiusOfCurvature, radialReserveX, radialReserveY, ...
        xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, ...
        xGridSimulationVec, yGridSimulationVec, zGridSimulationVec);
    if ~proceed
        disp('Simulation stopped by user.');
        return;
    end
end

% Define the medium parameters at the source position.
c0_source = MaterialMatrixSimulation.c0(1,1,1);
rho0_source = MaterialMatrixSimulation.rho0(1,1,1);

% Keep only the data needed after setup.
if keepOverlayData
    [MaterialMatrixForPlot, xGridMediumVecForPlot, yGridMediumVecForPlot, zGridMediumVecForPlot] = ...
        prepare_visualization_overlay_data(MaterialMatrix, xGridMediumVec, yGridMediumVec, zGridMediumVec, false);
else
    MaterialMatrixForPlot = [];
    xGridMediumVecForPlot = [];
    yGridMediumVecForPlot = [];
    zGridMediumVecForPlot = [];
end
clear MaterialMatrix xGridMediumVec yGridMediumVec zGridMediumVec

% Determine the type of the boundary condition and calculate the velocity hologram.
shiftHoloDist = shiftHoloDistInWl * c0_source / frequency;
if isSphericalSource
    transducerThickness = radiusOfCurvature - sqrt(radiusOfCurvature^2 - (aperture/2)^2);
    if izBoundaryCondition > izApex + round((transducerThickness + shiftHoloDist)/dz)
        VelocityHologram = fp_flat_boundary_condition_sf_vec_core(TransducerSf, soundSpeedContactMedium, densityContactMedium, xDDxCalculationFlag, ...
            xBoundaryCondition, yBoundaryCondition, zBoundaryCondition, dx, dy, ServiceParameters);
    else
        VelocityHologram = set_flat_boundary_condition_sf_vec_core(TransducerSf, soundSpeedContactMedium, densityContactMedium, xDDxCalculationFlag, ...
            xBoundaryCondition, yBoundaryCondition, zBoundaryCondition, dx, dy, shiftHoloDistInWl, ServiceParameters);
    end
else
    TransducerSf = resample_TransducerSf(TransducerSf, ...
        xBoundaryCondition, yBoundaryCondition, zBoundaryCondition, dx, dy);
    VelocityHologram = TransducerSf.complexVelocityAmplitude;
end

% =========================================================================
% RUN SIMULATION
% =========================================================================

% Define the medium structure for k-Wave.
medium = [];
medium.sound_speed = MaterialMatrixSimulation.c0;
medium.density = MaterialMatrixSimulation.rho0;
medium.alpha_coeff = MaterialMatrixSimulation.alpha;
medium.alpha_power = alphaPower;
clear MaterialMatrixSimulation

% Setup the time grid.
[kGridSimulation, dt, nTimeGrid, pointsPerPeriod] = setup_time_grid( ...
    kGridSimulation, medium, nxSimulationGrid, nySimulationGrid, nzSimulationGrid, dx, dy, dz, frequency, CFL);

% Define the source structure for k-Wave.
source = [];
source.p_mask = false(nxSimulationGrid, nySimulationGrid, nzSimulationGrid);
source.p_mask(:,:,1) = true;
source.p_amplitude = abs(VelocityHologram(:)) * c0_source * rho0_source;
source.p_phase = angle(VelocityHologram(:)); % rad
source.p_frequency = frequency;

% Define the sensor box from 1D grid vectors to avoid full 3D coordinate copies.
ixSensor = find((single(xFieldBegin) <= single(xGridSimulationVec)) & (single(xGridSimulationVec) <= single(xFieldEnd)));
iySensor = find((single(yFieldBegin) <= single(yGridSimulationVec)) & (single(yGridSimulationVec) <= single(yFieldEnd)));
izSensor = find((single(zFieldBegin) <= single(zGridSimulationVec)) & (single(zGridSimulationVec) <= single(zFieldEnd)));
if isempty(ixSensor) || isempty(iySensor) || isempty(izSensor)
    error('Sensor mask is empty: no grid points lie inside the field box [xFieldBegin,xFieldEnd] x [yFieldBegin,yFieldEnd] x [zFieldBegin,zFieldEnd]. Check field limits vs. simulation grid.');
end

nxSensor = numel(ixSensor);
nySensor = numel(iySensor);
nzSensor = numel(izSensor);

% Define the sensor structure for k-Wave.
sensor = [];
sensor.mask = false(nxSimulationGrid, nySimulationGrid, nzSimulationGrid);
sensor.mask(ixSensor, iySensor, izSensor) = true;
sensor.record = {'p'};
sensor.record_start_index = nTimeGrid - pointsPerPeriod + 1;

% Build sensor coordinates only for the requested output volume.
[xSensor, ySensor, zSensor] = ndgrid( ...
    xGridSimulationVec(ixSensor), ...
    yGridSimulationVec(iySensor), ...
    zGridSimulationVec(izSensor));
clear ixSensor iySensor izSensor

% Define additional input arguments for the simulation.
input_args = {'PMLSize', [xSizePML ySizePML zSizePML], 'PMLInside', false, ...
    'PMLAlpha', [xAlphaPML yAlphaPML zAlphaPML], 'BinaryPath', kWaveBinPath, 'BinaryName', kWaveBinName};

if strongMemorySavingMode
    input_args = [input_args, {'StrongSaveToDisk', true}];
end

if prepareSimulationOnly
    input_args = [input_args, {'DeleteData', false, 'PrepareCommandOnly', true}];
    if ~isempty(preparedSimulationDataPath)
        input_args = [input_args, {'DataPath', preparedSimulationDataPath}];
    end
elseif ~isempty(presimulatedOutputFile)
    input_args = [input_args, {'DeleteData', false, 'LoadOutputFile', presimulatedOutputFile}];
end

% Add output format if provided in the script inputs (C++/CUDA only).
if ~isempty(outputSingleFrequencyFormat)
    input_args = [input_args, {'OutputSingleFrequencyFormat', outputSingleFrequencyFormat}];
end

% Run the k-Wave simulation.
switch kWaveCalculationFlag
    case 'cpu'
        % C++
        sensorData = kspaceFirstOrder3DC(kGridSimulation, medium, source, sensor, input_args{:});
    case 'cuda'
        % C++/CUDA GPU
        sensorData = kspaceFirstOrder3DG(kGridSimulation, medium, source, sensor, input_args{:});
end
clear kGridSimulation medium source sensor input_args

if prepareSimulationOnly
    fprintf(['k-Wave input preparation is complete. The external solver was not started.\n' ...
        'Run the printed command from your system terminal to generate the output HDF5.\n' ...
        'Then set showPresimulatedData to the file below and rerun this script for plotting:\n' ...
        'simulatorInputs.showPresimulatedData = ''%s'';\n'], sensorData.output_filename);
    return;
end

% Extract the pressure field.
if isfield(sensorData, 'p') && size(sensorData.p, 2) == 1 && ~isreal(sensorData.p)
    % Output is already complex (amplitude * exp(1i*phase)) from single-frequency mode.
    complexPressure = sensorData.p;
else
    % Traditional mode: extract amplitude and phase from time series.
    [amp, phase] = extractAmpPhase(sensorData.p, 1/dt, frequency, 'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
    complexPressure = amp.*exp(1i*phase);
end
clear sensorData dt amp phase

complexPressure = reshape(complexPressure, nxSensor, nySensor, nzSensor);

simulatorData.VelocityHologram = VelocityHologram;
simulatorData.xGridSimulationVec = xGridSimulationVec;
simulatorData.yGridSimulationVec = yGridSimulationVec;
simulatorData.zGridSimulationVec = zGridSimulationVec;
simulatorData.xSensor = xSensor;
simulatorData.ySensor = ySensor;
simulatorData.zSensor = zSensor;
simulatorData.complexPressure = complexPressure;
simulatorData.MaterialMatrixForPlot = MaterialMatrixForPlot;
simulatorData.xGridMediumVecForPlot = xGridMediumVecForPlot;
simulatorData.yGridMediumVecForPlot = yGridMediumVecForPlot;
simulatorData.zGridMediumVecForPlot = zGridMediumVecForPlot;
simulatorData.xBoundaryCondition = xBoundaryCondition;
simulatorData.yBoundaryCondition = yBoundaryCondition;
simulatorData.zBoundaryCondition = zBoundaryCondition;
simulatorData.TransducerSf = TransducerSf;
simulatorData.frequency = frequency;
simulatorData.aperture = aperture;
simulatorData.apertureNominalX = apertureNominalX;
simulatorData.apertureNominalY = apertureNominalY;
simulatorData.radiusOfCurvature = radiusOfCurvature;
simulatorData.isSphericalSource = isSphericalSource;
simulatorData.dx = dx;
simulatorData.dz = dz;
end

function validate_TransducerSf_type(TransducerSf)
    if ~isfield(TransducerSf, 'type')
        return;
    end

    transducerType = TransducerSf.type;
    if isa(transducerType, 'string')
        transducerType = char(transducerType);
    end
    validTypes = {'standard_single_element', 'standard_array', 'custom'};
    if ~ischar(transducerType) || ~any(strcmp(transducerType, validTypes))
        error('TransducerSf.type must be ''standard_single_element'', ''standard_array'', or ''custom'' if present. Delete the field to use an untyped transducer.');
    end
end

function validate_water_test_transducer(TransducerSf)
    if ~isfield(TransducerSf, 'type')
        error('waterTest requires a standard single-element transducer. TransducerSf.type is missing.');
    end

    transducerType = TransducerSf.type;
    if isa(transducerType, 'string')
        transducerType = char(transducerType);
    end
    if ~ischar(transducerType) || ~strcmp(transducerType, 'standard_single_element')
        error('waterTest requires a standard single-element transducer. TransducerSf.type must be ''standard_single_element''.');
    end
end

function throwAutoSelectionError(causeME)
    ME = MException('xDDx:AutoSimulationDeviceFailed', ...
        ['Simulation failed while automatic device selection was enabled. ', ...
        'Check the underlying error for the actual cause. Set device options manually if the failure is related to device selection.']);
    ME = addCause(ME, causeME);
    throw(ME);
end

function tf = shouldSaveAutoSelection(autoDeviceSelection, simulatorData)
    tf = autoDeviceSelection.usesAuto && isfield(simulatorData, 'complexPressure');
end
