% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% ***Quick-Start Hologram Alignment and Back-Projection***
%
% Applicable for focused transducers shaped as a spherical cup and flat
% transducers.
%
% Quick-start demonstration of the toolbox’s capabilities for hologram
% back-projection with or without alignment in the single frequency or transient
% regime. It uses the minimum possible number of input parameters, while
% the majority of the simulation parameters are defined automatically.
% This can be useful for quick checks of the holography measurements. For
% precise results in terms of reconstructing the vibrational pattern of the
% transducer surface, use the step by step approach presented in the
% manual.
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
%     'zPosition': approximate z-position in m of the hologram, assuming
%                  zero at the apex of the transducer
%     'dx': x-step of the hologram Cartesian grid in m
%     'dy': y-step of the hologram Cartesian grid in m
%     'time' : vector with the time samples in s at each time grid nodes
%              for the recorded waveforms
%     'pressureWaveforms': 2D or 3D matrix with the acoustic pressure
%                          waveforms in Pa at corresponding grid nodes
%                          of the hologram
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
libraryDir = '..\..\xDDx_lib'; %lib directory

%I. SIMULATION REGIME
simulationDevice = 'cuda'; %simulation device: 'cuda' or 'cpu'
doTransientCase = false; %true for transient case, false for single-frequency case
doAlignment = false; %set true to perform alignment, or false to perform back-projection with no alignment

%II. INPUT TRANSIENT HOLOGRAM DATA
inputParametersFileName = '..\data_for_examples\flat\holo_data_flat_transient.mat'; %input data with transient hologram and medium parameters


%III. TRANSDUCER PARAMETERS
frequencyNominal = 1e6; %desired frequency of the output single-frequency hologram in Hz
apertureMin = 38e-3; % characteristic minimum aperture of the transducer in m
apertureMax = 38e-3; % characteristic maximum aperture of the transducer in m (apertureMax = apertureMin for circular transducers)
radiusOfCurvature = []; % (optional) nominal radius of curvature of the transducer in m. Omit this variable or define it as an empty value [] if your transducer is flat


%IV. 3D FOCAL LOBE PARAMETERS
focalLobeLevel = 0.6; %isolevel related to the pressure maximum for extracting the isosurface of the main focal lobe

% First iteration: coarse grid step
ppwRadialCoarse3D = 3; %number of points per wavelength (ppw) in the radial direction for identifying the main lobe with a coarse grid step
ppwAxialCoarse3D = 2; %number of points per wavelength (ppw) in the axial direction for identifying the main lobe with a coarse grid step

% Second iteration: fine grid step
ppwRadialFine3D = 20; %number of points per wavelength (ppw) in the radial direction for identifying the main lobe with a fine grid step
ppwAxialFine3D = 10; %number of points per wavelength (ppw) in the axial direction for identifying the main lobe with a fine grid step


%TECHNICAL PARAMETERS
errorPositioningWl = 10; %expected error of the positioning of the hologram in the axial and radial directions in wavelengths of frequencyNominal
minHoloShiftInWl = 4; %minimum shift in wavelengths for the rotated hologram
powerSignificanceLevel = 0.01; %minimum significant power level of a spectral component relative to the maximum power
levelSignalToSave = 0.1;  %minimum velocity waveform level relative to the maximum value to be considered as the beginning of the signal
showRingingDown = false;   %true to show all output time points, false to show only the part of the signal greater than levelSignalToSave
apertureReserve = 0.1; %fraction of aperture used to increase the nominal output window
saveTransducer = false; %true to save the boundary condition for the transducer at frequencyNominal as a 'TransducerSf' structure

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
focalLobeTransparency = 0.4; %transparency of the isosurface of the main focal lobe
figuresCascadeShift = 50; % cascade shift for each new figure in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
isSphericalSource = exist('radiusOfCurvature', 'var') ~= 0;
if isSphericalSource
    isSphericalSource = isSphericalSource*(~isempty(radiusOfCurvature));
end

addpath(genpath(libraryDir));

%Calculate supplementary parameters
radiusMax = 0.5*(1+apertureReserve)*max(apertureMin, apertureMax);

%Load input data
load(inputParametersFileName);

if isvector(HologramTr.xGrid)
    [HologramTr] = reshape_transducer_dim(HologramTr);
end

%Calculate the holo-shift and the positioning error in m
minHoloShift = minHoloShiftInWl*Medium.soundSpeed/frequencyNominal;
errorPositioning = errorPositioningWl*Medium.soundSpeed/frequencyNominal;

%Create input Geometry structure
Geometry = [];
if isSphericalSource
    Geometry.radiusOfCurvature = radiusOfCurvature;
end
Geometry.apertureMax = apertureMax;
Geometry.apertureMin = apertureMin;

%Get the single-sided spectrum of the input transient signal
expSign = 1;
[fArray, pSpectrum] = get_single_sided_spectrum(HologramTr.time, HologramTr.pressureWaveforms);

%Extract the nominal frequency frequencyNominal from the spectrum
[~, ifNominal] = min(abs(fArray - frequencyNominal));
frequency = fArray(ifNominal);

%Create input structure HologramSf for the signal component at
%the nominal frequency frequencyNominal
HologramSf = [];
HologramSf.expSign = expSign;
HologramSf.frequency = frequency;
HologramSf.xGrid = HologramTr.xGrid;
HologramSf.yGrid = HologramTr.yGrid;
HologramSf.zPosition = HologramTr.zPosition;
HologramSf.dx = HologramTr.dx;
HologramSf.dy = HologramTr.dy;
HologramSf.complexPressureAmplitude = squeeze((pSpectrum(:,:, ifNominal)));
HologramSf.errorPositioning = errorPositioning;

%Get default figure position
positionCurrent = get(groot,'DefaultFigurePosition');


[doAlignment, warningFlatTransducerIssue] = check_flat_transducer_issue(isSphericalSource, doAlignment);

if doAlignment

    %Find alignment parameters for the mispositioned hologram
    [xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, RotationLine, fvFocalLobeCoarse, fvFocalLobeFine] = auto_alignment(simulationDevice, Geometry,...
        HologramSf, Medium, ServiceParameters, ...
        focalLobeLevel, ppwRadialCoarse3D, ppwAxialCoarse3D, ppwRadialFine3D, ppwAxialFine3D, minHoloShift);

    %Plot the coarse-grid focal lobe isosurface
    fvFocalLobeCoarseToPlot = fvFocalLobeCoarse;
    fvFocalLobeCoarseToPlot.vertices = fvFocalLobeCoarseToPlot.vertices*1e3;
    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    patch(fvFocalLobeCoarseToPlot, 'FaceAlpha', focalLobeTransparency, 'EdgeColor','none', 'FaceColor', 'r');
    axis equal
    axis tight;
    camlight;
    view(3);
    title('Focal lobe coarse grid step');
    xlabel('{\itx}, mm');
    ylabel('{\ity}, mm');
    zlabel('{\itz}, mm');


    %Plot the fine-grid focal lobe isosurface, field maximums for each z-layer, and an OLS optimized symmetry line
    fvFocalLobeFineToPlot = fvFocalLobeFine;
    fvFocalLobeFineToPlot.vertices = fvFocalLobeFineToPlot.vertices*1e3;
    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    patch(fvFocalLobeFineToPlot, 'FaceAlpha', focalLobeTransparency, 'EdgeColor','none', 'FaceColor', [34 177 76]/256);
    hold on
    plot3(RotationLine.xyzInputPoints(1,:)*1e3,RotationLine.xyzInputPoints(2,:)*1e3,RotationLine.xyzInputPoints(3,:)*1e3,'o');
    plot3(RotationLine.xyzLine(1,:)*1e3,RotationLine.xyzLine(2,:)*1e3,RotationLine.xyzLine(3,:)*1e3,'r');
    axis equal;
    axis tight;
    view(3);
    camlight;
    title(['Angle with the mechanical axis is ' num2str(RotationLine.angleZ) '^{\circ}']);
    xlabel('{\itx}, mm');
    ylabel('{\ity}, mm');
    zlabel('{\itz}, mm');

    %Forward-project and rotate hologram
    SourceParameters.xGrid = HologramSf.xGrid;
    SourceParameters.yGrid = HologramSf.yGrid;
    SourceParameters.zGrid = zeros(size(HologramSf.xGrid));
    SourceParameters.dx = HologramSf.dx;
    SourceParameters.dy = HologramSf.dy;    
    SourceParameters.input = HologramSf.complexPressureAmplitude;
    regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points
    isTransient = false;
    [xHoloRotated, yHoloRotated, zHoloRotated, holoRotated] = rotate_hologram(xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, Medium, ServiceParameters, radiusOfCurvature, minHoloShift);

    xHoloToSimulate = xHoloRotated;
    yHoloToSimulate = yHoloRotated;
    zHoloToSimulate = zHoloRotated;
    complexPressureAmplitudeToSimulate = holoRotated;
else

    xHoloToSimulate = HologramSf.xGrid;
    yHoloToSimulate = HologramSf.yGrid;
    zHoloToSimulate = HologramSf.zPosition*ones(size(xHoloToSimulate));
    complexPressureAmplitudeToSimulate = HologramSf.complexPressureAmplitude;

end


%Calculate the default steps of the source grid
dxSource = Medium.soundSpeed/frequency/2;
dySource = Medium.soundSpeed/frequency/2;

if ~isOctave
    dxSource = round(dxSource,1,'significant');
    dySource = round(dxSource,1,'significant');
end

nxSource = round((1+apertureReserve)*apertureMin/dxSource);
nySource = round((1+apertureReserve)*apertureMax/dySource);

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

if ~doTransientCase
    %Back-project the rotated hologram to the surface of the source
    SourceParameters.xGrid = xSourceInCirc;
    SourceParameters.yGrid = ySourceInCirc;
    SourceParameters.zGrid = zSourceInCirc;
    SourceParameters.input = [];

    FieldParameters.xGrid = xHoloToSimulate;
    FieldParameters.yGrid = yHoloToSimulate;
    FieldParameters.zGrid = zHoloToSimulate;
    FieldParameters.input = complexPressureAmplitudeToSimulate;
    FieldParameters.dx = HologramTr.dx;
    FieldParameters.dy = HologramTr.dy;

    isTransient = false;

    if isSphericalSource
        regime = 2;  % 2 Back-projection: P on a plane --> V on a sphere
        [ vSourceInCirc ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);
    else
        regime = 1;  % 1 Back-projection: P on a plane --> V on a plane
        [ vSourceInCirc ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);
    end

    %Reshape the output velocity back at the rectangular grid for subsequent visualization
    vSource = zeros(size(xSource));
    vSource(xSource.^2 + ySource.^2 <= radiusMax^2) = vSourceInCirc;

    %Plot complex velocity amplitude and phase at the surface of the transducer
    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,abs(vSource)*1e3);
    set(gca,'YDir','normal');
    colormap('jet');
    axis equal;
    axis tight;
    title(['Vibrational velocity amplitude [mm/s] at the array surface @ ' num2str(frequency*1e-6) ' MHz']);
    xlabel('{\itx}, mm');
    ylabel('{\ity}, mm');
    if isOctave
        xlim(xSource([1 end])*1e3)
        ylim(ySource([1 end])*1e3)
    end
    colorbar;

    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,angle(vSource));
    set(gca,'YDir','normal');
    colormap('hsv');
    axis equal;
    axis tight;
    title(['Vibrational velocity phase [rad.] at the array surface @ ' num2str(frequency*1e-6) ' MHz']);
    xlabel('{\itx}, mm');
    ylabel('{\ity}, mm');
    if isOctave
        xlim(xSource([1 end])*1e3)
        ylim(ySource([1 end])*1e3)
    end
    colorbar;

else

    %Find the number of significant frequency samples for the transient
    %hologram
    timeStep = HologramTr.time(2) - HologramTr.time(1);
    numberOfTimePoints = size(HologramTr.pressureWaveforms,3);
    frequencyStep = (timeStep*numberOfTimePoints)^-1;

    maxFrequency = Medium.soundSpeed/2/max(HologramTr.dx, HologramTr.dy);
    maxNumFrequencySamples = round(maxFrequency/frequencyStep);

    maxFrequencySample = find_highest_significant_freq(frequencyStep, pSpectrum, powerSignificanceLevel, maxNumFrequencySamples, maxFrequency);
    minFrequencySample = find_lowest_significant_freq(pSpectrum, powerSignificanceLevel);

    if doAlignment
        %Calculate global coordinates of the rotated hologram
        [xHoloRotatedGlobal, yHoloRotatedGlobal, zHoloRotatedGlobal, zPositionHoloRotated] = calculate_rotated_coordinates(xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, HologramTr.xGrid, HologramTr.yGrid, radiusOfCurvature, minHoloShift);

        SourceParameters = [];
        SourceParameters.xGrid = HologramTr.xGrid;
        SourceParameters.yGrid = HologramTr.yGrid;
        SourceParameters.zGrid = zeros(size(HologramTr.xGrid));
        SourceParameters.dx = HologramTr.dx;
        SourceParameters.dy = HologramTr.dy;

        %Forward-project and align the transient hologram
        SourceParameters.input = pSpectrum(:,:,2:(maxFrequencySample+1));

        FieldParameters = [];
        FieldParameters.xGrid = xHoloRotatedGlobal;
        FieldParameters.yGrid = yHoloRotatedGlobal;
        FieldParameters.zGrid = zHoloRotatedGlobal;
        regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points

        isTransient = true; %set to true for the transient regime, false for the single-frequency regime

        [ pSpectrumRotated ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

        pSpectrumToSimulate = pSpectrumRotated;

    else

        pSpectrumToSimulate = pSpectrum(:,:,2:(maxFrequencySample+1));

        if minFrequencySample > 1
            pSpectrumToSimulate(:,:,1:(minFrequencySample-1)) = 0;
        elseif minFrequencySample <= 0
            error('minFrequencySample must be greater than zero!');
        end

    end

    %Back-project the rotated transient hologram to the surface of the source
    SourceParameters = [];
    SourceParameters.xGrid = xSourceInCirc;
    SourceParameters.yGrid = ySourceInCirc;
    SourceParameters.zGrid = zSourceInCirc;
    SourceParameters.input = [];

    FieldParameters = [];
    FieldParameters.xGrid = xHoloToSimulate;
    FieldParameters.yGrid = yHoloToSimulate;
    FieldParameters.zGrid = zHoloToSimulate;
    FieldParameters.input = pSpectrumToSimulate;
    FieldParameters.dx = HologramTr.dx;
    FieldParameters.dy = HologramTr.dy;

    isTransient = true; %set to true for the transient regime, false for the single-frequency regime

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
    vWaveformsSurf = real(ifft(vSpectrumSurfTotal,[],3));

    %Shift the signal to find beginning and end points basing on the levelSignalToSave value
    [vWaveformsSurf, iBegin, iEnd] = circshift_waveform(vWaveformsSurf,levelSignalToSave,showRingingDown);

    %GUI part to show results a frame-by-frame format and save them
    if ~isSphericalSource
        radiusOfCurvature = [];
    end

    gui_transient_frame_by_frame(xSource, ySource, vSpectrumSurfTransient, frequencyStep, minFrequencySample, size(vSpectrumSurfTransient,3), false, radiusOfCurvature, dxSource, dySource, expSign);
    gui_transient_frame_by_frame(xSource, ySource, vWaveformsSurf, timeStep, iBegin, iEnd, true, radiusOfCurvature, dxSource, dySource, expSign, figuresCascadeShift);

end

if ~isempty(warningFlatTransducerIssue)
    warning(warningFlatTransducerIssue);
end

if saveTransducer

    if doTransientCase
        iFrequency = round(frequencyNominal/frequencyStep);

        if(iFrequency > size(vSpectrumSurfTransient,3))
            iFrequency = size(vSpectrumSurfTransient,3);
        end

        if(iFrequency < 1)
            iFrequency = 1;
        end

        frequency = iFrequency*frequencyStep;
        vSource = squeeze(vSpectrumSurfTransient(:,:,iFrequency));
    end


    save_transducer(isSphericalSource, expSign, HologramSf.frequency, xSource, ySource, zSource, radiusOfCurvature, dxSource, dySource, vSource);

end

disp('xDDx simulation completed!');