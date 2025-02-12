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

regimesAllVector = 1:6;

if ismac
    error(errorMessages.noMac);
elseif isunix
    error(errorMessages.noLinux);
elseif ~ispc
    disp(errorMessages.noOS);
end


spatialTolerance = eps('single');

ServiceParameters.percentStep = 1; %default size of a percent step for the percent counter

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

libFolder = cd(tempdir);

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

cd(cppExeFolder);
cppExeName = ['rayleigh_' simulationPostfix '.exe'];
status = system(cppExeName);
if status ~= 0
    error(errorMessages.exe);
end

cd(tempdir);

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

delete(BinFileNames.xSource);
delete(BinFileNames.ySource);
delete(BinFileNames.zSource);
delete(BinFileNames.xField);
delete(BinFileNames.yField);
delete(BinFileNames.zField);
delete(BinFileNames.vectorParam);
delete(BinFileNames.reInput);
delete(BinFileNames.imInput);
delete(BinFileNames.reOutput);
delete(BinFileNames.imOutput);

cd(libFolder);

end

