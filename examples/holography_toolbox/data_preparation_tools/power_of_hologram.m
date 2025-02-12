% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%This tool is designed to calculate and plot a power curve for different 
% spectral components of a transient hologram. It helps determine the
% maximum significant frequency for transient hologram post-processing. 
% The tool uses an angular spectrum-based method for power calculation
% in a single-frequency planar scan 
%
% INPUT DATA FORMAT FOR A TRANSIENT HOLOGRAM: a MAT file
% 'inputParametersFileName' with the following variables
%
% 'Medium' struct with fields:
%     'soundSpeed': sound speed in m/s
%     'density': density in kg/m^3
%
% 'HologramTr' struct with fields:
%     'xGrid': vector or matrix with x-coordinates of the Cartesian grid in
%              m at each grid node of the hologram
%     'yGrid': vector or matrix with y-coordinates of the Cartesian grid in
%              m at each grid node of the hologram
%     'dx': x-step of the hologram Cartesian grid in m
%     'dy': y-step of the hologram Cartesian grid in m
%     'zPosition': approximate z-position in m of the hologram, assuming 
%                  zero at the apex of the transducer
%     'time' : vector with the time samples in s at each time grid nodes 
%              for the recorded waveforms
%     'pressureWaveforms': 3D matrix with the acoustic pressure 
%                          waveforms in Pa at corresponding grid nodes of the hologram
%                          size(pressureWaveforms) = [size(xGrid, 1) size(xGrid, 2) length(time)]
% Here, x and y are the transverse coordinates of the transducer, and z is 
% the coordinate along the direction of the beam. See details in the manual.


clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%EDITABLE CODE%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
libraryDir = '..\..\..\xDDx_lib'; %lib directory

%I. INPUT TRANSIENT HOLOGRAM DATA
inputParametersFileName = '..\..\data_for_examples\spherical\holo_data_spherical_transient.mat'; %input data with transient hologram and medium parameters
stepDim1 = 0.4e-3; %step size in m along the 1st dimension of the HologramTr.pressureWaveforms 3D matrix
stepDim2 = 0.4e-3; %step size in m along the 2nd dimension of the HologramTr.pressureWaveforms 3D matrix

%SERVICE PARAMETERS
figuresCascadeShift = 50; % cascade shift for each new figure in pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(libraryDir));

%Load input data
load(inputParametersFileName);

if isvector(HologramTr.xGrid)
  error('This tool works only for 2D hologram grids!');
end

%Perform fft to extract spectral components for signals at all hologram points
timeStep = HologramTr.time(2) - HologramTr.time(1);
numberOfTimePoints = size(HologramTr.pressureWaveforms,3);
frequencyStep = (timeStep*numberOfTimePoints)^-1;
pSpectrum = fft(HologramTr.pressureWaveforms,[],3);
pSpectrum = pSpectrum/numberOfTimePoints;
if mod(numberOfTimePoints,2) == 0
    pSpectrum = pSpectrum(:,:,1:(numberOfTimePoints/2+1));
    pSpectrum(:,:,2:(end-1)) = 2*pSpectrum(:,:,2:(end-1));
else
    pSpectrum = pSpectrum(:,:,1:((numberOfTimePoints+1)/2));
    pSpectrum(:,:,2:end) = 2*pSpectrum(:,:,2:end);
end
pSpectrum = pSpectrum(:,:,2:end);

frequencyArray = (1:size(pSpectrum,3))*frequencyStep;

%Calculate the power
powerAtFrequency = zeros(size(frequencyArray));
for iFreq = 1:length(frequencyArray)
    powerAtFrequency(iFreq) = holo_pressure_to_power(pSpectrum(:,:,iFreq), stepDim1, stepDim2, frequencyArray(iFreq), Medium);
end

%Show results
positionCurrent = get(groot,'DefaultFigurePosition');

%Show results
figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
plot(frequencyArray*1e-6, powerAtFrequency);
xlabel('frequency, MHz');
ylabel('power of the spectral component, W');

figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
plot(powerAtFrequency);
xlabel('frequency sample');
ylabel('power of the spectral component, W');

disp('xDDx simulation completed!');