% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% ***Transient Hologram Rotation and Back-Projection***
%
% This tool rotates and shifts a measured transient hologram so that the 
% aligned hologram is perpendicular to the transducer's axis of symmetry, 
% then back-projects it to the transducer's surface. The script reads 
% predetermined position parameters from the file 
% 'alignmentParametersFileName'. This file can either be saved using the 
% Automatic Alignment Tool (bp_auto_alignment_sf) or constructed manually.
% The output frame-by-frame GUI displays results in both frequency and 
% time domains, with an option to render the results in video format.
%
% Note: transient simulations may require more processing time compared to 
% single-frequency simulations.
%
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
%     'pressureWaveforms': 2D or 3D matrix with the acoustic pressure
%                          waveforms in Pa at corresponding grid nodes of the hologram
%                          if xGrid/yGrid is a 2D matrix, then
%                          size(pressureWaveforms) = [size(xGrid, 1) size(xGrid, 2) length(time)]
%                          if xGrid/yGrid is a vector, then
%                          size(pressureWaveforms) = [length(xGrid) length(time)]
%
% POSITION PARAMETERS: a MAT file 'alignmentParametersFileName' with the
% following variables
%     'zMaxAcoustCoord': the acoustical z-coordinate of the maximum 
%        pressure point. For strongly focused transducers, this coordinate 
%        equals the radius of curvature.
%     'xMaxMechCoord', 'yMaxMechCoord', 'zMaxMechCoord': coordinates of 
%        the pressure maximum point in the coordinates of the measured 
%        hologram with the origin at its center ("mechanical coordinates")
%     'directionVectorMechCoord': unit vector in the "mechanical 
%        coordinates" directed along the acoustical axis of the transducer 
%        towards the focus point
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
simulationDevice = 'cuda'; %simulation device: 'cuda' or 'cpu'

%II. INPUT PARAMETERS
inputParametersFileName = '..\..\data_for_examples\spherical\holo_data_spherical_transient.mat'; %input data with transient hologram, and medium parameters
frequencyCentral = 1.25e6; %approximate central frequency of the transducer in Hz

%III. ALIGNMENT PARAMETERS
alignmentParametersFileName = '..\..\data_for_examples\spherical\alig_param_spherical.mat'; %data structure with position parameters (can be saved using 'run_auto_alignment_sf.m')

%IV. SOURCE PARAMETERS
nxSource = 251; %number of points for the Cartesian grid of the source along the x-dimension
nySource = 251; %number of points for the Cartesian grid of the source along the y-dimension
dxSource = 0.4e-3; %x grid step in m
dySource = 0.4e-3; %y grid step in m

%V. FREQUENCY DOMAIN PARAMETERS
powerSignificanceLevel = 0.01; %minimum significant power level of a spectral component relative to the maximum power (used in 'auto' regime)
minFrequencySample = 'auto'; %frequency samples ranging (minFrequencySample:maxFrequencySample)*frequencyStep, where frequencyStep = 1 / "time window for HologramTr.time"
maxFrequencySample = 'auto'; %set 'auto' to determine the values automatically based on the powerSignificanceLevel

%VI. TECHNICAL PARAMETERS
minHoloShiftInWl = 4; %minimum shift in wavelengths of the central frequency for the rotated hologram
levelSignalToSave = 0.1; %minimum velocity waveform level relative to the maximum velocity to be saved in the output video
showRingingDown = false; %true to show all output time points, false to show only the part of the signal greater than levelSignalToSave
apertureReserve = 0.2; %fraction of aperture used to increase the nominal output window

%VII. OUTPUT VIDEOS PARAMETERS
saveSurfSpectrum = false;  %set true to save a video with the back-projected velocity amplitudes at different frequencies
outputVideoSurfSpectrum = 'out_v_spec_source.mp4'; %output filename
frameRateVideoSurfSpectrum = 4; %frame rate in the output video

saveSurfSignal = false; %set true to save a video with the back-projected velocity waveforms at the surface of the transducer
outputVideoSurfSignal = 'out_v_sgnl_source.mp4'; %output filename
downsampleCoeffTime = 4; %downsampling factor for the output back-projected velocity waveform related to the input 'HologramTr.time' samples
frameRateVideoSurfSignal = 10; %frame rate in the output video

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; % cascade shift for each new figure in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(libraryDir));

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
[saveSurfSpectrum, saveSurfSignal] = check_octave_based_video_issues(isOctave, saveSurfSpectrum, saveSurfSignal);

%Load input data
load(inputParametersFileName);
load(alignmentParametersFileName);

if isvector(HologramTr.xGrid)
    [HologramTr] = reshape_transducer_dim(HologramTr);
end

%Calculate supplementary parameters
radiusMax = 0.5*(1+apertureReserve)*max(Geometry.apertureMin, Geometry.apertureMax);

%Set geometric parameters
if isfield(Geometry, 'radiusOfCurvature')
    if ~isempty(Geometry.radiusOfCurvature)
        radiusOfCurvature = zMaxAcoustCoord;
    end
end
isSphericalSource = exist('radiusOfCurvature', 'var') ~= 0;
if isSphericalSource
  isSphericalSource = isSphericalSource*(~isempty(radiusOfCurvature));
end

%Perform fft to extract spectral components for signals at all hologram points
expSign = 1;
timeStep = HologramTr.time(2) - HologramTr.time(1);
numberOfTimePoints = size(HologramTr.pressureWaveforms,3);
frequencyStep = (timeStep*numberOfTimePoints)^-1;
pSpectrum = fft(HologramTr.pressureWaveforms,[],3);

%Set 'auto' frequency range if needed
if strcmpi(minFrequencySample,'auto')
   [minFrequencySample] = find_lowest_significant_freq(pSpectrum, powerSignificanceLevel);
end

if strcmpi(maxFrequencySample,'auto')
   maxFrequency = Medium.soundSpeed/2/max(HologramTr.dx, HologramTr.dy);
   maxNumFrequencySamples = round(maxFrequency/frequencyStep);
   [maxFrequencySample] = find_highest_significant_freq(frequencyStep, pSpectrum, powerSignificanceLevel, maxNumFrequencySamples, maxFrequency);
end

%Renormalize “MATLAB/Octave fft” into “single-sided fft”
pSpectrum = pSpectrum/numberOfTimePoints;
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

%Calculate the grid for the rotated hologram using preloaded ALIGNMENT PARAMETERS
xHolo = HologramTr.xGrid;
yHolo = HologramTr.yGrid;

%Calculate the holo-shift in m
minHoloShift = minHoloShiftInWl*Medium.soundSpeed/frequencyCentral; %minimum shift in m for the rotated hologram
[xHoloRotatedGlobal, yHoloRotatedGlobal, zHoloRotatedGlobal, zPositionHoloRotated] = calculate_rotated_coordinates(xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, xHolo, yHolo, zMaxAcoustCoord, minHoloShift);

SourceParameters.xGrid = xHolo;
SourceParameters.yGrid = yHolo;
SourceParameters.zGrid = zeros(size(xHolo));
SourceParameters.dx = HologramTr.dx;
SourceParameters.dy = HologramTr.dy;

%Forward-project and align transient hologram
SourceParameters.input = pSpectrum;

FieldParameters.xGrid = xHoloRotatedGlobal;
FieldParameters.yGrid = yHoloRotatedGlobal;
FieldParameters.zGrid = zHoloRotatedGlobal;
regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points

isTransient = true; %set to true for the transient regime, false for the single-frequency regime

[ pSpectrumRotated ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

xHoloRotated = xHolo;
yHoloRotated = yHolo;
zHoloRotated = zPositionHoloRotated*ones(size(xHoloRotated));

%Generate a Cartesian grid on a spherical surface of the transducer
[xSource, ySource] = meshgrid((-(nxSource-1)/2:(nxSource-1)/2)*dxSource,...
    (-(nySource-1)/2:(nySource-1)/2)*dySource);

if isSphericalSource
    zSource = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xSource.^2 - ySource.^2);
else
    zSource = zeros(size(xSource));
end

%Select the source grid nodes inside circle of interest (radiusMax)
xSourceInCirc = xSource(xSource.^2 + ySource.^2 <= radiusMax^2);
ySourceInCirc = ySource(xSource.^2 + ySource.^2 <= radiusMax^2);
zSourceInCirc = zSource(xSource.^2 + ySource.^2 <= radiusMax^2);

%Back-project the rotated transient hologram to the surface of the source
SourceParameters.xGrid = xSourceInCirc;
SourceParameters.yGrid = ySourceInCirc;
SourceParameters.zGrid = zSourceInCirc;
SourceParameters.input = [];

FieldParameters.xGrid = xHoloRotated;
FieldParameters.yGrid = yHoloRotated;
FieldParameters.zGrid = zHoloRotated;
FieldParameters.input = pSpectrumRotated;
FieldParameters.dx = HologramTr.dx;
FieldParameters.dy = HologramTr.dy;

if isSphericalSource
    regime = 2; % 2 Back-projection: P on a plane --> V on a sphere
   [ vSpectrumSurfTransientInCirc ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);
else
    regime = 1; % 1 Back-projection: P on a plane --> V on a plane
   [ vSpectrumSurfTransientInCirc ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);
end

%Reshape the output velocity back at the rectangular grid for subsequent visualization
vSpectrumSurfTransient = zeros(size(xSource,1)*size(xSource,2), size(vSpectrumSurfTransientInCirc,3));
vSpectrumSurfTransient(xSource.^2 + ySource.^2 <= radiusMax^2,:) = vSpectrumSurfTransientInCirc;
vSpectrumSurfTransient = reshape(vSpectrumSurfTransient,[size(xSource) size(vSpectrumSurfTransientInCirc,3)]);

%Renormalize single-sided fft into MATLAB/Octave fft
vSpectrumSurfTransientMatlabCalibration = vSpectrumSurfTransient*numberOfTimePoints;
if mod(numberOfTimePoints,2) == 0
    vSpectrumSurfTransientMatlabCalibration(:,:,2:(end-1)) = vSpectrumSurfTransientMatlabCalibration(:,:,2:(end-1))/2;
else
    vSpectrumSurfTransientMatlabCalibration(:,:,2:end) = vSpectrumSurfTransientMatlabCalibration(:,:,2:end)/2;
end

%Zero-pad the obtained spectrum to maintain the same time step of the output as in input time grid 'HologramTr.time'
vSpectrumSurfTotal = cat(3, zeros(size(vSpectrumSurfTransientMatlabCalibration,1),size(vSpectrumSurfTransientMatlabCalibration,2)),...
    vSpectrumSurfTransientMatlabCalibration,...
    zeros(size(vSpectrumSurfTransientMatlabCalibration,1),size(vSpectrumSurfTransientMatlabCalibration,2),(numberOfTimePoints-2*maxFrequencySample-1)),...
    flip(conj(vSpectrumSurfTransientMatlabCalibration),3));

%Return back in time domain
vWaveformsSurf = real(ifft(vSpectrumSurfTotal,[],3));

%Save a video with single-frequency holograms in the range of frequencies (1:(numberOfFrequencySamples-1))*frequencyStep
if saveSurfSpectrum
    figure
    v = VideoWriter(outputVideoSurfSpectrum, 'MPEG-4');
    set(v,'FrameRate',frameRateVideoSurfSpectrum);
    set(v,'Quality',100);
    open(v);
    for i = minFrequencySample:1:size(vSpectrumSurfTransient,3)
        imagesc(xSource([1 end])*1e3, ySource([1 end])*1e3 , squeeze(abs(vSpectrumSurfTransient(:,:,i)))*1e3);
        set(gca,'YDir','normal');

        axis equal

        colormap jet
        h = colorbar;
        axis tight;

        title(['Vibrational velocity amplitude [mm/s] at {\itf} =' num2str(i*frequencyStep*1e-6, '%1.4f') ' MHz']);
        xlabel('x, mm');
        ylabel('y, mm');
        frame_out = getframe(gcf);
        writeVideo(v,frame_out);

    end
    close(v);
end

%Shift the signal to find beginning and end points basing on the levelSignalToSave value
[vWaveformsSurf, iBegin, iEnd] = circshift_waveform(vWaveformsSurf,levelSignalToSave,showRingingDown);

%Save a video of the vibrational velocity at the surface of the transducer at different time moments
if saveSurfSignal

    vLimits = [min(real(vWaveformsSurf(:))) max(vWaveformsSurf(:))]*1e3;

    figure
    v = VideoWriter(outputVideoSurfSignal, 'MPEG-4');
    set(v,'FrameRate',frameRateVideoSurfSignal);
    set(v,'Quality',100);
    open(v);
    for i = iBegin:downsampleCoeffTime:iEnd
        imagesc(xSource([1 end])*1e3, ySource([1 end])*1e3 , squeeze(vWaveformsSurf(:,:,i))*1e3);
        set(gca,'YDir','normal');

        axis equal

        caxis(vLimits);

        colormap jet
        h = colorbar;
        axis tight;

        title(['Vibrational velocity [mm/s] at {\itt} =' num2str((i-iBegin)*timeStep*1e6, '%1.4f') ' us']);
        xlabel('x, mm');
        ylabel('y, mm');
        frame_out = getframe(gcf);
        writeVideo(v,frame_out);

    end
    close(v);

end

%GUI part to show results a frame-by-frame format and save them
if ~isSphericalSource
  radiusOfCurvature = [];
end

gui_transient_frame_by_frame(xSource, ySource, vSpectrumSurfTransient, frequencyStep, minFrequencySample, size(vSpectrumSurfTransient,3), false, radiusOfCurvature, dxSource, dySource, expSign);
gui_transient_frame_by_frame(xSource, ySource, vWaveformsSurf, timeStep, iBegin, iEnd, true, radiusOfCurvature, dxSource, dySource,  expSign, figuresCascadeShift);

disp('xDDx simulation completed!');