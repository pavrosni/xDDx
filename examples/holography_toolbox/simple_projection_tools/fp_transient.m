% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% ***Transient Hologram Forward-Projection***
%
% This tool is designed for the forward projection of a transient hologram 
% to obtain the pressure waveform at user-defined points.
%
% Note: transient simulations may require more processing time compared to
% single-frequency simulations.
%
% INPUT DATA FORMAT FOR A TRANSIENT HOLOGRAM: a MAT file
% 'inputParametersFileName' with the following variables
%
% 'Medium' struct with fields:
%     'soundSpeed': sound speed in m/s
%     'density': density in kg/m^3
%
% 'HologramTr' struct with fields:
%     'xGrid': vector or matrix with the x-coordinates in m at each grid
%     node of the hologram  (Cartesian grid only)
%     'yGrid': vector or matrix with the y-coordinates in m at
%     corresponding grid nodes of the hologram  (Cartesian grid only)
%     'zPosition' (optional field): z-position in m of the hologram,
%     assuming zero at the apex of the transducer  (Cartesian grid only)
%     'dx': x-step of the hologram Cartesian grid in m
%     'dy': y-step of the hologram Cartesian grid in m
%     'time' : vector with the time samples in s at each time grid nodes
%     for the recorded waveforms 
%     'pressureWaveforms': 2D or 3D matrix with the acoustic pressure
%     waveforms in Pa at corresponding grid nodes of the hologram 
%                          if xGrid/yGrid is a 2D matrix, then
%                          size(pressureWaveforms) = [size(xGrid, 1) size(xGrid, 2) length(time)] 
%                          if xGrid/yGrid is a vector, then
%                          size(pressureWaveforms) = [length(xGrid) length(time)] 
%
% Here, x and y are the transverse coordinates of the transducer, and z is 
% the coordinate along the direction of the beam. See details in the manual.

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%EDITABLE CODE%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
libraryDir = '..\..\..\xDDx_lib'; %lib directory

%I. SIMULATION REGIME
simulationDevice = 'cpu'; %simulation device: 'cuda' or 'cpu'

%II. INPUT PARAMETERS
inputParametersFileName = '..\..\data_for_examples\spherical\holo_data_spherical_transient.mat'; %input data with transducer, hologram, and medium parameters

%III. POINTS FOR FIELD SIMULATION
xField = [0 -0.25e-3 1e-3]; %x, y, and z coordinates in m of the field points of interest
yField = [0 0 1e-3];
zField = [88e-3 88.5e-3 90e-3];

%V. FREQUENCY DOMAIN PARAMETERS
powerSignificanceLevel = 0.01; %minimum significant power level of a spectral component relative to the maximum power
minFrequencySample = 'auto'; %frequency samples ranging (minFrequencySample:maxFrequencySample)*frequencyStep, where frequencyStep = 1 / "time window for HologramTr.time"
maxFrequencySample = 'auto'; %set 'auto' to determine the values automatically based on the powerSignificanceLevel

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; % cascade shift for each new figure in pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(libraryDir));

%Load input data
load(inputParametersFileName);

if isvector(HologramTr.xGrid)
    [HologramTr] = reshape_transducer_dim(HologramTr);
end

%Perform fft to extract spectral components for signals at all hologram points
expSign = 1;
timeStep = HologramTr.time(2) - HologramTr.time(1);
numberOfTimePoints = size(HologramTr.pressureWaveforms,3);
frequencyStep = (timeStep*numberOfTimePoints)^-1;
pSpectrum = fft(HologramTr.pressureWaveforms,[],3);
pSpectrum = pSpectrum/numberOfTimePoints;

%Set 'auto' frequency range if needed
if strcmpi(minFrequencySample,'auto')
   [minFrequencySample] = find_lowest_significant_freq(pSpectrum, powerSignificanceLevel);
end

if strcmpi(maxFrequencySample,'auto')
   maxFrequency = Medium.soundSpeed/2/max(HologramTr.dx, HologramTr.dy);
   maxNumFrequencySamples = round(maxFrequency/frequencyStep);
   [maxFrequencySample] = find_highest_significant_freq(frequencyStep, pSpectrum, powerSignificanceLevel, maxNumFrequencySamples, maxFrequency);
end

if mod(numberOfTimePoints,2) == 0
    pSpectrum = pSpectrum(:,:,1:(numberOfTimePoints/2+1));
    pSpectrum(:,:,2:(end-1)) = 2*pSpectrum(:,:,2:(end-1));
else
    pSpectrum = pSpectrum(:,:,1:((numberOfTimePoints+1)/2));
    pSpectrum(:,:,2:end) = 2*pSpectrum(:,:,2:end);
end

pSpectrum = pSpectrum(:,:,2:(maxFrequencySample+1));

if minFrequencySample > 1
    pSpectrum(:,:,1:(minFrequencySample-1)) = 0;
elseif minFrequencySample <= 0
    error('minFrequencySample must be greater than zero!');
end



%Forward-project the transient hologram
isTransient = true; 

xHolo = HologramTr.xGrid;
yHolo = HologramTr.yGrid;

SourceParameters.xGrid = xHolo;
SourceParameters.yGrid = yHolo;
SourceParameters.zGrid = HologramTr.zPosition*ones(size(xHolo));
SourceParameters.dx = HologramTr.dx;
SourceParameters.dy = HologramTr.dy;
SourceParameters.input = pSpectrum;

FieldParameters.xGrid = xField;
FieldParameters.yGrid = yField;
FieldParameters.zGrid = zField;
regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points
[ pSpectrumField ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);


%Renormalize single-sided fft into MATLAB/Octave fft
pSpectrumFieldMatlabCalibration = pSpectrumField*numberOfTimePoints;
if mod(numberOfTimePoints,2) == 0
    pSpectrumFieldMatlabCalibration(:,:,2:(end-1)) = pSpectrumFieldMatlabCalibration(:,:,2:(end-1))/2;
else
    pSpectrumFieldMatlabCalibration(:,:,2:end) = pSpectrumFieldMatlabCalibration(:,:,2:end)/2;
end

%Zero-pad the obtained spectrum to maintain the same time step of the output as in input time grid 'HologramTr.time'
pSpectrumFieldTotal = cat(3, zeros(size(pSpectrumFieldMatlabCalibration,1),size(pSpectrumFieldMatlabCalibration,2)),...
    pSpectrumFieldMatlabCalibration,...
    zeros(size(pSpectrumFieldMatlabCalibration,1),size(pSpectrumFieldMatlabCalibration,2),(numberOfTimePoints-2*maxFrequencySample-1)),...
    flip(conj(pSpectrumFieldMatlabCalibration),3));

%Return back in time domain
pFieldWaveforms = squeeze(real(ifft(pSpectrumFieldTotal,[],3)));

%Plot the result

if isvector(pFieldWaveforms)
    pFieldWaveforms = reshape(pFieldWaveforms, [1 length(pFieldWaveforms)]);
end

positionCurrent = get(groot,'DefaultFigurePosition');
for iField = 1:numel(xField)

    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];

    plot(HologramTr.time*1e6, pFieldWaveforms(iField,:)*1e-3, 'k');
    xlim([HologramTr.time(1), HologramTr.time(end)]*1e6);
    xlabel('Time, {\mu}s');
    ylabel('Pressure, kPa');
    title(['Pressure waveform  (x = ' num2str(xField(iField)*1e3) ' mm, ' 'y = ' num2str(yField(iField)*1e3) ' mm, ' 'z = ' num2str(zField(iField)*1e3) ' mm)' ]);
    grid on;
    grid minor;
end

disp('xDDx simulation completed!');

