% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% 
% USAGE:
% [ outputField ] = rayleigh_simulator(expSign, frequencyParameter, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature)
% 
% 
% INPUTS:
%     expSign               +1 or -1, depending on the exponent sign convention exp(+ 1i * omega * t) or exp(- 1i * omega * t).
% 	                        E.g., the "fft" function in MATLAB utilizes the  exp(+ 1i * omega * t) convention, so in this case, expSign is +1
% 
%     frequencyParameter    frequencyParameter = frequency of the transducer in Hz for a single-frequency hologram, or frequencyParameter = frequencyStep for a transient hologram
% 
%     regime                simulation regime number from 1 to 6, see details below*
% 
%     simulationDevice      'cuda' - perform simulation using CUDA compatible videocard (GPU)
%                           'cpu'  - perform simulation using the central processor (CPU) of the computer
%                           'auto' - automatically select 'cuda' when available, otherwise 'cpu'
% 
%     isTransient           true for the transient regime, false for the single-frequency regime
% 
%     SourceParameters      input or output (depending on the regime) structure that describes the Source, see details below**
% 
%     FieldParameters       input or output (depending on the regime) structure that describes the Field, see details below**
% 
%     Medium                structure with medium parameters, see details below***
% 
%     radiusOfCurvature     (optional, can be omitted if unnecessary) radius of curvature in m of the input spherical surface of integration in regimes 3 and 6. Omit this variable if your transducer is flat. 
%
%     ServiceParameters     (optional, can be omitted if unnecessary) structure with service parameters, see details below. If omitted, the default ServiceParameters are set.****
% 
% 
% OUTPUT:
%     outputField           matrix of acoustic pressure/vibrational velocity complex amplitude at the Source/Field surface (depending on the regime), see details regarding the matrix size below*****  
% 
% 
% *REGIMES:
% 1 Back-projection: P on a plane --> V on a plane
% 2 Back-projection: P on a plane --> V on a sphere
% 3 Forward-projection: V on a sphere --> P at an arbitrary set of points
% 4 Forward-projection: V on a plane --> P at an arbitrary set of points
% 5 Forward-projection: P on a plane --> P at an arbitrary set of points
% 6 Back-projection: V on a sphere --> P on a planes 
% 
% 
% **SIMULATION PARAMETERS FORMAT:
% 'SourceParameters' and 'FieldParameters' structures that may contain the following fields
%     'xGrid' (necessary field): vector or matrix with the x-coordinates in m at each grid node of the Source or Field region (Cartesian grid only)
%     'yGrid' (necessary field): vector or matrix with the y-coordinates in m at each grid node of the Source or Field region (Cartesian grid only)
%     'zGrid' (necessary field): vector or matrix with the z-coordinates in m at each grid node of the Source or Field region (Cartesian grid only)
%     'dx' (optional field): x-step of the Source or Field Cartesian grid in m
%     'dy' (optional field): y-step of the Source or Field Cartesian grid in m
%     'input'(optional field):  input complex pressure or velocity amplitude area for integration surface in Pa or m/s, see details regarding the matrix size below*****
% 
% 
% ***MEDIUM PARAMETERS:
% 'Medium' struct with fields:
%     'soundSpeed': sound speed in m/s
%     'density': density in kg/m^3
%
% 
% ****SERVICE PARAMETERS:
% 'ServiceParameters'  (optional) struct with fields:
%     'threadsPerBlockGPU': number of threads per block for GPU if applicable. Default value is 128.
%     'cpuArchitecture': CPU executable/image architecture ('auto', 'avx512', 'avx2', 'avx', 'sse2', or 'arm64'). Default value is 'auto'.
%     'cudaVersion': CUDA executable/image version ('auto', 'cuda11', 'cuda12', 11, or 12). Default value is 'auto'.
%
% 
% *****INPUT/OUTPUT FIELD MATRIX:
% for a single-frequency hologram (isTransient = false):
% input/output complex amplitudes are given for each node of the the input/output grid 
% i.e. size(input) = size(xGrid)  
% 
% for a transient hologram (isTransient = true):
% input/output complex amplitudes are given for each node of the the input/output grid in a range of frequencies (1:numberOfFrequencySamples)*frequencyStep
% i.e. size(input) = [size(xGrid) numberOfFrequencySamples] 

function [ outputField ] = rayleigh_simulator(expSign, frequencyParameter, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, varargin)

[BinFileNames] = load_bin_file_names;
[errorMessages] = load_error_messages;
[originalSimulationDevice, autoDeviceSelection] = resolve_auto_simulation_device(simulationDevice);
simulationDevice = originalSimulationDevice;

if (length(varargin) > 2) || isempty(varargin)
     error(errorMessages.simInput);
end

for iArg = 1:length(varargin)
  if isstruct(varargin{iArg})
      ServiceParameters = varargin{iArg};
  elseif isnumeric(varargin{iArg})
      radiusOfCurvature = varargin{iArg};
  else
      error(errorMessages.simInput);
  end
end

if exist('ServiceParameters', 'var') == 0
    ServiceParameters = [];
    ServiceParameters.threadsPerBlockGPU = 128; 
end    

if exist('radiusOfCurvature', 'var') == 0
    radiusOfCurvature = 1;
end    

[ServiceParameters, usedCachedCudaVersion] = apply_cached_cuda_version(ServiceParameters, autoDeviceSelection, simulationDevice);
selectedCudaVersion = '';

regimesAllVector = 1:6;

% On Linux/Mac, simulations run via Docker. On Windows, local executables are used
% unless XDDX_USE_DOCKER_ON_WINDOWS is enabled.
if ismac && strcmp(simulationDevice, 'cuda')
    error('CUDA simulation is not supported on macOS. Use CPU mode instead (simulationDevice = ''cpu'').');
end
if ~ispc && ~ismac && ~isunix
    disp(errorMessages.noOS);
end

privateWorkDir = '';
try

spatialTolerance = eps('single');

ServiceParameters.percentStep = 5; %default size of a percent step for the percent counter

[inputRayleigh, surfElementArea, isInputSource, isInputField] = check_input_errors(simulationDevice, isTransient, regime, regimesAllVector, errorMessages,SourceParameters,FieldParameters, spatialTolerance, radiusOfCurvature);

if isTransient && isvector(SourceParameters.xGrid)
    if isInputSource
        SourceParameters.xGrid = reshape(SourceParameters.xGrid,[1 numel(SourceParameters.xGrid)]);
        SourceParameters.yGrid = reshape(SourceParameters.yGrid,[1 numel(SourceParameters.yGrid)]);
        SourceParameters.zGrid = reshape(SourceParameters.zGrid,[1 numel(SourceParameters.zGrid)]);
        
        squeezedInputRayleigh = squeeze(inputRayleigh);
        inputRayleigh = reshape(inputRayleigh, [1 size(squeezedInputRayleigh,1) size(squeezedInputRayleigh,2)]);
    end
end

if isTransient && isvector(FieldParameters.xGrid)
    if isInputField
        FieldParameters.xGrid = reshape(FieldParameters.xGrid,[1 numel(FieldParameters.xGrid)]);
        FieldParameters.yGrid = reshape(FieldParameters.yGrid,[1 numel(FieldParameters.yGrid)]);
        FieldParameters.zGrid = reshape(FieldParameters.zGrid,[1 numel(FieldParameters.zGrid)]);
        
        squeezedInputRayleigh = squeeze(inputRayleigh);
        inputRayleigh = reshape(inputRayleigh, [1 size(squeezedInputRayleigh,1) size(squeezedInputRayleigh,2)]);
    end
end

vecInputParam = [Medium.density Medium.soundSpeed frequencyParameter expSign regime radiusOfCurvature surfElementArea ServiceParameters.percentStep];

if strcmp(simulationDevice,'cuda')
    vecInputParam = [vecInputParam ServiceParameters.threadsPerBlockGPU];
end

simulationPostfix = simulationDevice;
if isTransient
    simulationPostfix = [simulationPostfix '_transient'];   
end


cppExeFolder = fullfile(fileparts(mfilename('fullpath')), 'rayleigh_cpp', simulationPostfix);

rayleighWorkDir = tempdir;
if ismac
    % Avoid rapidly reusing shared file paths between successive Docker runs.
    rayleighWorkDir = tempname(tempdir);
    [created, message] = mkdir(rayleighWorkDir);
    if ~created
        error('xDDx:TempDirectory', 'Cannot create %s: %s', rayleighWorkDir, message);
    end
    privateWorkDir = rayleighWorkDir;
end
libFolder = cd(rayleighWorkDir);

write_matrix_bin(BinFileNames.vectorParam, vecInputParam);



write_matrix_bin(BinFileNames.xSource, SourceParameters.xGrid);
write_matrix_bin(BinFileNames.ySource, SourceParameters.yGrid);
write_matrix_bin(BinFileNames.zSource, SourceParameters.zGrid);

if ~ismatrix(FieldParameters.xGrid)
    write_matrix_bin(BinFileNames.xField, FieldParameters.xGrid(:));
    write_matrix_bin(BinFileNames.yField, FieldParameters.yGrid(:));
    write_matrix_bin(BinFileNames.zField, FieldParameters.zGrid(:));
else
    write_matrix_bin(BinFileNames.xField, FieldParameters.xGrid);
    write_matrix_bin(BinFileNames.yField, FieldParameters.yGrid);
    write_matrix_bin(BinFileNames.zField, FieldParameters.zGrid);
end
  
write_matrix_bin(BinFileNames.reInput, real(inputRayleigh));
write_matrix_bin(BinFileNames.imInput, imag(inputRayleigh));

useDocker = (ismac || isunix) || force_docker_win();
if useDocker
    % Pre-flight: ensure Docker is installed and usable (e.g. no sudo required on Linux)
    dockerCheck = check_docker_ready();
    % On macOS, MATLAB/system() may run with a PATH that doesn't include
    % Homebrew/Docker Desktop locations (common cause: `zsh: command not found: docker`).
    % If Docker CLI isn't found in PATH, try to locate it automatically and
    % update PATH for subsequent docker commands.
    % Use strfind/lower for older MATLAB compatibility (avoid `contains`).
    if ~dockerCheck.ok && ~isempty(strfind(lower(dockerCheck.message), 'not found in path'))
        dockerPathFix = ensure_docker_cli_in_path();
        if dockerPathFix.didChangePath
            dockerCheck = check_docker_ready();
        end
    end
    if ~dockerCheck.ok
        nl = char(10);
        error(['Docker check failed: ' dockerCheck.message nl nl 'Hint: ' dockerCheck.hint]);
    end
    % Run simulation via Docker (Linux/Mac, or Windows when explicitly requested)
    dockerConfig = xddx_docker_config();
    imageName = get_xddx_docker_image(simulationDevice, ServiceParameters, isTransient);
    selectedCudaVersion = get_cuda_version_from_image_name(imageName);
    fullImage = [dockerConfig.dockerUsername '/' imageName];
    runImage = ensure_docker_image_ready(fullImage, imageName, dockerConfig);
    tempDir = rayleighWorkDir;
    if tempDir(end) == filesep
        tempDir = tempDir(1:end-1);
    end
    volumeArg = ['-v "' tempDir ':/Temp"'];
    runRm = '--rm';
    if strcmp(simulationDevice, 'cuda')
        dockerCmd = sprintf('docker run %s --gpus all %s %s', runRm, volumeArg, runImage);
    else
        dockerCmd = sprintf('docker run %s %s %s', runRm, volumeArg, runImage);
    end
    [status, cmdout] = run_command_live(dockerCmd);
    if status ~= 0
        errMsg = sprintf(errorMessages.dockerRun, runImage);
        runHint = get_docker_run_failure_hint(cmdout);
        if ~isempty(runHint)
            errMsg = [errMsg char(10) char(10) 'Hint: ' runHint];
        end
        error(errMsg);
    end
else
    % Run via local exe (Windows)
    cd(cppExeFolder);
    [cppExeName, executableVariant] = get_local_cpp_exe_name(simulationPostfix, ServiceParameters);
    if startsWith(executableVariant, 'cuda')
        selectedCudaVersion = executableVariant;
    end
    cppExePath = fullfile(cppExeFolder, cppExeName);
    if exist(cppExePath, 'file') ~= 2
        error(['Rayleigh Integral Simulator executable is missing: ' cppExePath]);
    end
    status = system(quote_shell_arg(cppExePath));
    if status ~= 0
        error(get_local_cpp_exe_failure_message(errorMessages.exe, executableVariant));
    end
end

cd(rayleighWorkDir);

testReOutput = fopen(BinFileNames.reOutput);
testImOutput = fopen(BinFileNames.imOutput);

if (testReOutput == -1)||(testImOutput == -1)
    error(errorMessages.noCuda);
end

fclose(testReOutput);
fclose(testImOutput);

clear testFid;

outputField = read_matrix_bin(BinFileNames.reOutput) + 1i*read_matrix_bin(BinFileNames.imOutput);

if isempty(outputField) || any(isnan(outputField(:))) || (max(abs(outputField(:))) < eps)
    error(errorMessages.noCuda);
end

if ~ismatrix(FieldParameters.xGrid)
    if ~isTransient
    outputField = reshape(outputField, size(FieldParameters.xGrid));
    else
    outputField = reshape(outputField, [size(FieldParameters.xGrid) size(outputField, 3)]);
    end
end

% Mac files stay together until success so failures retain diagnostic inputs.
if ~ismac
    % Clean up temp files (Docker may remove or chown some on Linux).
    delete_if_exists(BinFileNames.xSource);
    delete_if_exists(BinFileNames.ySource);
    delete_if_exists(BinFileNames.zSource);
    delete_if_exists(BinFileNames.xField);
    delete_if_exists(BinFileNames.yField);
    delete_if_exists(BinFileNames.zField);
    delete_if_exists(BinFileNames.vectorParam);
    % Linux Docker files can be root-owned; overwrite these on the next run.
    if ~(isunix && useDocker)
        delete_if_exists(BinFileNames.reInput);
        delete_if_exists(BinFileNames.imInput);
        delete_if_exists(BinFileNames.reOutput);
        delete_if_exists(BinFileNames.imOutput);
    end
end

cd(libFolder);

if should_save_auto_selection(autoDeviceSelection, outputField)
    autoDeviceSelection = set_successful_auto_cuda_version(autoDeviceSelection, ...
        simulationDevice, selectedCudaVersion, usedCachedCudaVersion);
    resolve_auto_simulation_device(simulationDevice, ...
        'AutoDeviceSelection', autoDeviceSelection, ...
        'SaveSelection', true);
end

if ~isempty(privateWorkDir)
    % Remove only this invocation's directory after leaving it and reading outputs.
    try
        [removed, cleanupMessage] = rmdir(privateWorkDir, 's');
    catch cleanupError
        removed = false;
        cleanupMessage = cleanupError.message;
    end
    if ~removed
        warning('xDDx:TempCleanup', ...
            'Simulation completed, but temporary files remain in %s: %s', ...
            privateWorkDir, cleanupMessage);
    end
end

catch ME
    if exist('libFolder', 'var') == 1
        try
            cd(libFolder);
        catch
        end
    end
    if ~isempty(privateWorkDir) && exist(privateWorkDir, 'dir') == 7
        fprintf(2, 'Rayleigh diagnostic files retained in: %s\n', privateWorkDir);
    end
    outputField = handle_auto_simulation_failure(ME, autoDeviceSelection, ...
        expSign, frequencyParameter, regime, isTransient, ...
        SourceParameters, FieldParameters, Medium, varargin{:});
end

end

function outputField = handle_auto_simulation_failure(causeME, autoDeviceSelection, ...
    expSign, frequencyParameter, regime, isTransient, SourceParameters, FieldParameters, Medium, varargin)

if isempty(autoDeviceSelection) || ~isfield(autoDeviceSelection, 'usesAuto') || ~autoDeviceSelection.usesAuto
    rethrow(causeME);
end

if autoDeviceSelection.usedCache
    resolve_auto_simulation_device('auto', ...
        'AutoDeviceSelection', autoDeviceSelection, ...
        'ClearCache', true);
    [retrySimulationDevice, retryAutoDeviceSelection] = resolve_auto_simulation_device('auto', ...
        'ForceDetect', true);
    try
        outputField = rayleigh_simulator(expSign, frequencyParameter, regime, ...
            retrySimulationDevice, isTransient, SourceParameters, FieldParameters, Medium, varargin{:});
        if should_save_auto_selection(retryAutoDeviceSelection, outputField)
            resolve_auto_simulation_device(retrySimulationDevice, ...
                'AutoDeviceSelection', retryAutoDeviceSelection, ...
                'SaveSelection', true);
        end
        return;
    catch retryME
        throw_auto_selection_error(retryME);
    end
end

throw_auto_selection_error(causeME);
end

function throw_auto_selection_error(causeME)
ME = MException('xDDx:RayleighAutoSimulationDeviceFailed', ...
    ['Simulation failed while automatic device selection was enabled. ', ...
    'Check the underlying error for the actual cause. Set device options manually if the failure is related to device selection.']);
ME = addCause(ME, causeME);
throw(ME);
end

function tf = should_save_auto_selection(autoDeviceSelection, outputField)
tf = ~isempty(autoDeviceSelection) ...
    && isfield(autoDeviceSelection, 'usesAuto') ...
    && autoDeviceSelection.usesAuto ...
    && ~isempty(outputField);
end

function on = force_docker_win()
% Use Docker on Windows when the canonical or legacy override is enabled.
on = use_xddx_docker_on_windows();
end

function [ServiceParameters, usedCachedCudaVersion] = apply_cached_cuda_version(ServiceParameters, autoDeviceSelection, simulationDevice)
usedCachedCudaVersion = false;

if ~strcmp(simulationDevice, 'cuda') ...
        || isempty(autoDeviceSelection) ...
        || ~isfield(autoDeviceSelection, 'usedCache') ...
        || ~autoDeviceSelection.usedCache ...
        || ~isfield(autoDeviceSelection, 'selectedCudaVersion') ...
        || isempty(autoDeviceSelection.selectedCudaVersion) ...
        || has_explicit_cuda_version(ServiceParameters)
    return;
end

cachedCudaVersion = char(autoDeviceSelection.selectedCudaVersion);
if ~is_valid_cuda_version(cachedCudaVersion)
    return;
end

if ~isstruct(ServiceParameters)
    ServiceParameters = struct();
end

ServiceParameters.cudaVersion = cachedCudaVersion;
usedCachedCudaVersion = true;
end

function tf = has_explicit_cuda_version(ServiceParameters)
tf = false;
if ~isstruct(ServiceParameters) || ~isfield(ServiceParameters, 'cudaVersion')
    return;
end

requestedCudaVersion = ServiceParameters.cudaVersion;
if isnumeric(requestedCudaVersion)
    tf = true;
    return;
elseif isa(requestedCudaVersion, 'string')
    requestedCudaVersion = char(requestedCudaVersion);
end

if ~ischar(requestedCudaVersion)
    tf = true;
    return;
end

tf = ~strcmpi(strtrim(requestedCudaVersion), 'auto');
end

function tf = is_valid_cuda_version(cudaVersion)
tf = any(strcmp(char(cudaVersion), {'cuda11', 'cuda12'}));
end

function cudaVersion = get_cuda_version_from_image_name(imageName)
cudaVersion = '';
tokens = regexp(imageName, '(cuda11|cuda12)$', 'tokens', 'once');
if ~isempty(tokens)
    cudaVersion = tokens{1};
end
end

function autoDeviceSelection = set_successful_auto_cuda_version(autoDeviceSelection, simulationDevice, selectedCudaVersion, usedCachedCudaVersion)
if isempty(autoDeviceSelection) || ~isfield(autoDeviceSelection, 'usesAuto') || ~autoDeviceSelection.usesAuto
    return;
end

if strcmp(simulationDevice, 'cuda') && is_valid_cuda_version(selectedCudaVersion)
    autoDeviceSelection.selectedCudaVersion = selectedCudaVersion;
elseif isfield(autoDeviceSelection, 'selectedCudaVersion') && ~usedCachedCudaVersion
    autoDeviceSelection.selectedCudaVersion = '';
end
end

function [cppExeName, executableVariant] = get_local_cpp_exe_name(simulationPostfix, ServiceParameters)
baseExeName = ['rayleigh_' simulationPostfix];
executableVariant = '';

if startsWith(simulationPostfix, 'cuda')
    [cudaVersion, hasExplicitCudaVersion] = get_xddx_cuda_version(ServiceParameters);
    candidateCudaVersions = get_local_cuda_version_candidates(cudaVersion, hasExplicitCudaVersion);

    for iCudaVersion = 1:numel(candidateCudaVersions)
        candidateExeName = [baseExeName '-' candidateCudaVersions{iCudaVersion} '.exe'];
        if exist(fullfile(pwd, candidateExeName), 'file') == 2
            cppExeName = candidateExeName;
            executableVariant = candidateCudaVersions{iCudaVersion};
            return;
        end
    end

    legacyExeName = [baseExeName '.exe'];
    if exist(fullfile(pwd, legacyExeName), 'file') == 2
        cppExeName = legacyExeName;
        executableVariant = '';
        return;
    end

    cppExeName = [baseExeName '-' cudaVersion '.exe'];
    executableVariant = cudaVersion;
    return;
end

if ~startsWith(simulationPostfix, 'cpu')
    cppExeName = [baseExeName '.exe'];
    return;
end

[executableVariant, hasExplicitArchitecture] = get_service_cpu_architecture(ServiceParameters);
candidateArchitectures = get_local_cpu_architecture_candidates(executableVariant, hasExplicitArchitecture);

for iArchitecture = 1:numel(candidateArchitectures)
    candidateExeName = [baseExeName '-' candidateArchitectures{iArchitecture} '.exe'];
    if exist(fullfile(pwd, candidateExeName), 'file') == 2
        cppExeName = candidateExeName;
        executableVariant = candidateArchitectures{iArchitecture};
        return;
    end
end

legacyExeName = [baseExeName '.exe'];
if exist(fullfile(pwd, legacyExeName), 'file') == 2
    cppExeName = legacyExeName;
    executableVariant = '';
    return;
end

cppExeName = [baseExeName '-' executableVariant '.exe'];
end

function candidateCudaVersions = get_local_cuda_version_candidates(cudaVersion, hasExplicitCudaVersion)
if hasExplicitCudaVersion
    candidateCudaVersions = {cudaVersion};
elseif strcmp(cudaVersion, 'cuda12')
    candidateCudaVersions = {'cuda12', 'cuda11'};
else
    candidateCudaVersions = {'cuda11', 'cuda12'};
end
end

function [cpuArchitecture, hasExplicitArchitecture] = get_service_cpu_architecture(ServiceParameters)
[cpuArchitecture, hasExplicitArchitecture] = get_xddx_cpu_architecture(ServiceParameters);
end

function candidateArchitectures = get_local_cpu_architecture_candidates(cpuArchitecture, hasExplicitArchitecture)
if strcmp(cpuArchitecture, 'arm64')
    candidateArchitectures = {'arm64'};
    return;
end

if ispc && ~hasExplicitArchitecture
    candidateArchitectures = {'avx2', 'avx', 'sse2', 'avx512'};
    return;
end

architectureOrder = {'avx512', 'avx2', 'avx', 'sse2'};
startIdx = find(strcmp(architectureOrder, cpuArchitecture), 1);
if isempty(startIdx)
    candidateArchitectures = {cpuArchitecture};
else
    candidateArchitectures = architectureOrder(startIdx:end);
end
end

function errorMessage = get_local_cpp_exe_failure_message(baseErrorMessage, cpuArchitecture)
errorMessage = baseErrorMessage;
if isempty(cpuArchitecture)
    return;
end

if strcmp(cpuArchitecture, 'avx512')
    errorMessage = [baseErrorMessage char(10) char(10) ...
        'The AVX512 CPU executable failed. If this computer does not support AVX512, try ServiceParameters.cpuArchitecture = ''avx2''.'];
elseif strcmp(cpuArchitecture, 'avx2')
    errorMessage = [baseErrorMessage char(10) char(10) ...
        'The AVX2 CPU executable failed. If this computer does not support AVX2, try ServiceParameters.cpuArchitecture = ''avx''.'];
elseif strcmp(cpuArchitecture, 'avx')
    errorMessage = [baseErrorMessage char(10) char(10) ...
        'The AVX CPU executable failed. If this computer does not support AVX, try ServiceParameters.cpuArchitecture = ''sse2''.'];
elseif strcmp(cpuArchitecture, 'sse2')
    errorMessage = [baseErrorMessage char(10) char(10) ...
        'The SSE2 CPU executable failed. This computer may be unsupported or the executable may be missing dependencies.'];
elseif strcmp(cpuArchitecture, 'cuda12')
    errorMessage = [baseErrorMessage char(10) char(10) ...
        'The CUDA 12 executable failed. First try ServiceParameters.cudaVersion = ''cuda11''. If CUDA 11 also fails, run the simulation in CPU mode (simulationDevice = ''cpu'').'];
elseif strcmp(cpuArchitecture, 'cuda11')
    errorMessage = [baseErrorMessage char(10) char(10) ...
        'The CUDA 11 executable failed. Run the simulation in CPU mode (simulationDevice = ''cpu''). If this machine should support newer CUDA drivers, update the NVIDIA driver and try CUDA again.'];
end
end

function delete_if_exists(fname)
% Delete file if it exists; avoid "File not found or permission denied" warnings
% (e.g. when Docker has removed input files or left them root-owned on Linux).
if exist(fname, 'file')
    prev = warning('off', 'all');
    try
        delete(fname);
    catch
    end
    warning(prev);
end
end

function runImage = ensure_docker_image_ready(fullImage, imageName, dockerConfig)
% Use a local image immediately if present. Refresh it only when due, and
% fall back to the local cached image if pulling fails (e.g. offline use).

runImage = find_local_docker_image(fullImage, imageName);

% If the configured Docker Hub username is wrong but the desired image is
% already loaded locally under another repository name, use it immediately.
if ~isempty(runImage) && ~strcmp(runImage, fullImage)
    fprintf('Using local Docker image: %s\n', runImage);
    return;
end

localImageExists = ~isempty(runImage);
shouldPull = ~localImageExists || docker_pull_due(fullImage, dockerConfig.imageUpdatePeriodDays);

if ~shouldPull
    return;
end

fprintf('Checking Docker image: %s\n', fullImage);
[pullStatus, pullOutput] = run_command_live(sprintf('docker pull %s', fullImage));
record_docker_pull_attempt(fullImage, pullStatus == 0);

if pullStatus == 0
    runImage = fullImage;
    return;
end

runImage = find_local_docker_image(fullImage, imageName);
if ~isempty(runImage)
    fprintf('Docker pull failed, using local Docker image: %s\n', runImage);
    return;
end

nl = char(10);
errMsg = [ ...
    'Docker image is not available locally and could not be pulled:' nl ...
    fullImage nl nl ...
    'This usually means there is no internet connection, Docker Hub is unavailable, or the Docker username/image name is incorrect.' nl ...
    'Please connect to the internet and run "docker pull ' fullImage '", then rerun the simulation.' ...
    ];
pullHint = get_docker_run_failure_hint(pullOutput);
if ~isempty(pullHint)
    errMsg = [errMsg nl nl 'Hint: ' pullHint];
end
error(errMsg);
end

function localImage = find_local_docker_image(fullImage, imageName)
localImage = '';

if docker_image_exists(fullImage)
    localImage = fullImage;
    return;
end

if docker_image_exists(imageName)
    localImage = imageName;
    return;
end

[status, imageList] = system('docker images --format "{{.Repository}}:{{.Tag}}"');
if status ~= 0 || isempty(imageList)
    return;
end

lines = regexp(strtrim(imageList), '\r\n|\n|\r', 'split');
for iLine = 1:numel(lines)
    candidate = strtrim(lines{iLine});
    if isempty(candidate) || ~isempty(strfind(candidate, '<none>'))
        continue;
    end

    colonIdx = find(candidate == ':', 1, 'last');
    if isempty(colonIdx)
        repository = candidate;
    else
        repository = candidate(1:colonIdx-1);
    end

    if strcmp(repository, imageName) || string_ends_with(repository, ['/' imageName])
        localImage = candidate;
        return;
    end
end
end

function tf = docker_image_exists(fullImage)
inspectCommand = sprintf('docker image inspect %s %s', fullImage, get_shell_null_redirect());
status = system(inspectCommand);
tf = (status == 0);
end

function tf = string_ends_with(textValue, suffix)
if length(textValue) < length(suffix)
    tf = false;
    return;
end

tf = strcmp(textValue(end-length(suffix)+1:end), suffix);
end

function due = docker_pull_due(fullImage, updatePeriodDays)
state = load_docker_pull_state();
idx = find(strcmp({state.images.name}, fullImage), 1);

if isempty(idx)
    due = true;
    return;
end

lastAttempt = state.images(idx).lastAttempt;
if isempty(lastAttempt)
    due = true;
    return;
end

due = datetime('now', 'TimeZone', 'local') >= lastAttempt + days(updatePeriodDays);
end

function record_docker_pull_attempt(fullImage, wasSuccessful)
state = load_docker_pull_state();
idx = find(strcmp({state.images.name}, fullImage), 1);
nowTime = datetime('now', 'TimeZone', 'local');

if isempty(idx)
    idx = numel(state.images) + 1;
    state.images(idx).name = fullImage;
    state.images(idx).lastAttempt = [];
    state.images(idx).lastSuccess = [];
end

state.images(idx).lastAttempt = nowTime;
if wasSuccessful
    state.images(idx).lastSuccess = nowTime;
end

save_docker_pull_state(state);
end

function state = load_docker_pull_state()
stateFile = get_docker_pull_state_file();
state = struct('images', struct('name', {}, 'lastAttempt', {}, 'lastSuccess', {}));

if exist(stateFile, 'file') ~= 2
    return;
end

try
    loadedState = load(stateFile, 'state');
    if isfield(loadedState, 'state') && isfield(loadedState.state, 'images')
        state = loadedState.state;
    end
catch
end
end

function save_docker_pull_state(state)
stateFile = get_docker_pull_state_file();
try
    save(stateFile, 'state');
catch
end
end

function stateFile = get_docker_pull_state_file()
stateFile = fullfile(prefdir, 'xddx_docker_pull_state.mat');
end

function [status, cmdout] = run_command_live(command)
% Run a shell command while streaming its merged stdout/stderr to MATLAB.
% This preserves progress output in the Command Window and still returns
% the full text for post-run error analysis.

if ~usejava('jvm')
    [status, cmdout] = system(command);
    if ~isempty(cmdout)
        fprintf('%s', cmdout);
    end
    return;
end

shellCommand = get_shell_command(command);
commandArray = javaArray('java.lang.String', numel(shellCommand));
for iCmd = 1:numel(shellCommand)
    commandArray(iCmd) = java.lang.String(shellCommand{iCmd});
end

builder = java.lang.ProcessBuilder(commandArray);
builder.redirectErrorStream(true);

try
    process = builder.start();
catch
    [status, cmdout] = system(command);
    if ~isempty(cmdout)
        fprintf('%s', cmdout);
    end
    return;
end

reader = java.io.InputStreamReader(process.getInputStream());
capturedOutput = java.lang.StringBuilder();
charCount = 0;

while true
    nextChar = reader.read();
    if nextChar == -1
        break;
    end

    nextCharMatlab = char(nextChar);
    fprintf('%c', nextCharMatlab);
    capturedOutput.append(nextCharMatlab);

    charCount = charCount + 1;
    if mod(charCount, 256) == 0
        drawnow;
    end
end

reader.close();
status = process.waitFor();
drawnow;

cmdout = char(capturedOutput.toString());
end

function shellCommand = get_shell_command(command)
if ispc
    shellCommand = {'cmd.exe', '/c', command};
else
    shellCommand = {'/bin/sh', '-c', command};
end
end

function nullRedirect = get_shell_null_redirect()
if ispc
    nullRedirect = '>nul 2>&1';
else
    nullRedirect = '>/dev/null 2>&1';
end
end

function quotedArg = quote_shell_arg(arg)
quotedArg = ['"' strrep(arg, '"', '\"') '"'];
end

function fix = ensure_docker_cli_in_path()
% Try to auto-fix missing docker CLI on macOS/Linux by updating PATH.
% This is safe to run multiple times.
fix = struct('didChangePath', false, 'addedDirs', {{}});

% Request targets macOS, but this helper is harmless on other unix-like systems.
if ~(ismac || isunix)
    return;
end

% Candidate locations for the `docker` executable:
% - Homebrew (Apple Silicon): /opt/homebrew/bin/docker
% - Homebrew (Intel): /usr/local/bin/docker
% - Docker Desktop bundle:
%     /Applications/Docker.app/Contents/Resources/bin/docker
% - Some setups symlink docker into /usr/bin/docker
candidates = { ...
    '/opt/homebrew/bin/docker', ...
    '/usr/local/bin/docker', ...
    '/Applications/Docker.app/Contents/Resources/bin/docker', ...
    '/usr/bin/docker' ...
    };

oldPath = getenv('PATH');
if isempty(oldPath)
    oldPath = '';
end
% On macOS/Linux PATH separator is ":"; avoid `strsplit` for older MATLAB.
pathDirs = regexp(oldPath, ':', 'split');

for i = 1:numel(candidates)
    dockerExe = candidates{i};
    if exist(dockerExe, 'file') ~= 2
        continue;
    end

    dockerDir = fileparts(dockerExe);
    if any(strcmp(pathDirs, dockerDir))
        continue;
    end

    setenv('PATH', [dockerDir pathsep oldPath]);
    fix.didChangePath = true;
    fix.addedDirs{end+1} = dockerDir; %#ok<AGROW>

    % Prepend first working candidate; keep it simple and predictable.
    break;
end
end
