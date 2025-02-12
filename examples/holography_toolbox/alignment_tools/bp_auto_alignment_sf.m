% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% ***Automatic Single-Frequency Hologram Rotation and Backprojection***
%
% Applicable for focused transducers shaped as a spherical cup with a known
% nominal radius of curvature.
%
% Assumptions:
%  - The transducer's field has a symmetrical focal lobe to find the axis
%    of symmetry.
%  - The focal maximum is at the center of curvature (true for strongly
%    focused transducers).
%
% Automatically rotates and shifts a measured single-frequency hologram so
% that the rotated hologram is perpendicular to the axis of symmetry of the
% transducer and backprojects it to its surface.
%
% The output GUI figure can be used to change the radius of curvature of
% the transducer and re-backproject the field. The actual radius of
% curvature corresponds to the sharp backprojected image. Save the
% backprojected complex velocity amplitude field and the parameters
% describing the rotated hologram using the GUI.
%
% INPUT DATA FORMAT FOR A SINGLE-FREQUENCY HOLOGRAM: a MAT or XLSX file
% 'inputParametersFileName' with the following variables
%
% 'Geometry' struct with fields:
%   'radiusOfCurvature': nominal radius of curvature of the transducer in m
%   'apertureMin': characteristic minimum aperture of the transducer in m     
%   'apertureMax': characteristic maximum aperture of the transducer in m
%                  (apertureMax = apertureMin for circular transducers)     
%
% 'Medium' struct with fields:
%   'soundSpeed': sound speed in m/s
%   'density': density in kg/m^3
%
% 'HologramSf' struct with fields:
%   'expSign': +1 or -1, depending on the exponent sign convention
%       exp(+ 1i * omega * t) or exp(- 1i * omega * t). E.g., the "fft"
%       function in MATLAB utilizes the  exp(+ 1i * omega * t) convention,
%       so in this case, expSign is +1
%   'frequency': the operating frequency of the transducer in Hz
%   'xGrid': vector or matrix with the x-coordinates in m at each grid node
%       of the hologram (Cartesian grid only)
%   'yGrid': vector or matrix with the y-coordinates in m at corresponding
%       grid nodes of the hologram (Cartesian grid only)
%   'zPosition': z-position in m of the hologram, assuming zero at the apex
%       of the transducer (Cartesian grid only)
%   'dx': x-step of the hologram Cartesian grid in m
%   'dy': y-step of the hologram Cartesian grid in m
%   'complexPressureAmplitude': vector or matrix with the complex pressure
%       amplitude values in Pa at corresponding grid nodes of the hologram
%
%  XLSX template path is 
%  'xDDx\examples\data_for_examples\xlsx_templates\hologram_sf.xls'
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
inputParametersFileName = '..\..\data_for_examples\spherical\holo_data_spherical_sf.mat'; %input data with transducer, single-frequency hologram, and medium parameters in '.mat' or '.xlsx' format


%III. SOURCE PARAMETERS
nxSource = 201; %number of points for the Cartesian grid of the source along the x-dimension
nySource = 201; %number of points for the Cartesian grid of the source along the y-dimension
dxSource = 0.5e-3; %x grid step in m
dySource = 0.5e-3; %y grid step in m
xZeroPhase = -25e-3; %(optional parameter) x-coordinate of the zero-phase point on the surface of the source
yZeroPhase = -25e-3; %(optional parameter) y-coordinate of the zero-phase point on the surface of the source

%IV. 3D FOCAL LOBE PARAMETERS
focalLobeLevel = 0.6; %isolevel related to the pressure maximum for extracting the isosurface of the main focal lobe


%V. 3D FOCAL LOBE PARAMETERS
% First iteration: coarse grid step
ppwRadialCoarse3D = 3; %number of points per wavelength (ppw) in the radial direction for identifying the main lobe with a coarse grid step
ppwAxialCoarse3D = 2; %number of points per wavelength (ppw) in the axial direction for identifying the main lobe with a coarse grid step

% Second iteration: fine grid step
ppwRadialFine3D = 20; %number of points per wavelength (ppw) in the radial direction for identifying the main lobe with a fine grid step
ppwAxialFine3D = 10; %number of points per wavelength (ppw) in the axial direction for identifying the main lobe with a fine grid step

%VI. TECHNICAL PARAMETERS
errorPositioningWl = 10; %expected error of the positioning of the hologram in the axial and radial directions in wavelengths of HologramSf.frequency
minHoloShiftInWl = 4; %minimum shift of the central frequency for the rotated hologram in wavelengths of HologramSf.frequency
apertureReserve = 0.2; %fraction of aperture used to increase the nominal output window

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
focalLobeTransparency = 0.4; %transparency of the isosurface of the main focal lobe
figuresCascadeShift = 50; % cascade shift for each new figure in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

addpath(genpath(libraryDir));

%Load input data
[~, ~, inputExt] = fileparts(inputParametersFileName);
if strcmpi(inputExt, '.mat')
    load(inputParametersFileName);
else
    [Geometry, HologramSf, Medium] = read_hologram_sf_from_xls(inputParametersFileName);
end

%Calculate the holo-shift and the positioning error in m
frequency = HologramSf.frequency;
minHoloShift = minHoloShiftInWl*Medium.soundSpeed/frequency;
errorPositioning = errorPositioningWl*Medium.soundSpeed/frequency;

%Calculate and define supplementary parameters
radiusMax = 0.5*(1+apertureReserve)*max(Geometry.apertureMin, Geometry.apertureMax);

%Find alignment parameters for the mispositioned hologram
HologramSf.errorPositioning = errorPositioning;
[xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, RotationLine, fvFocalLobeCoarse, fvFocalLobeFine] = auto_alignment(simulationDevice, Geometry,...
    HologramSf, Medium, ServiceParameters, ...
    focalLobeLevel, ppwRadialCoarse3D, ppwAxialCoarse3D, ppwRadialFine3D, ppwAxialFine3D, minHoloShift);

%Plot the coarse-grid focal lobe isosurface
positionCurrent = get(groot,'DefaultFigurePosition');
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
expSign = HologramSf.expSign;

SourceParameters.xGrid = HologramSf.xGrid;
SourceParameters.yGrid = HologramSf.yGrid;
SourceParameters.zGrid = zeros(size(HologramSf.xGrid));
SourceParameters.dx = HologramSf.dx;
SourceParameters.dy = HologramSf.dy;
SourceParameters.input = HologramSf.complexPressureAmplitude;
regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points

isTransient = false;

radiusOfCurvature = Geometry.radiusOfCurvature;

[xHoloRotated, yHoloRotated, zHoloRotated, holoRotated] = rotate_hologram(xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, Medium, ServiceParameters, radiusOfCurvature, minHoloShift);

%Generate a Cartesian grid on a spherical surface of the transducer
[xSource, ySource] = meshgrid((-(nxSource-1)/2:(nxSource-1)/2)*dxSource,...
    (-(nySource-1)/2:(nySource-1)/2)*dySource);
zSource = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xSource.^2 - ySource.^2);

%Select the source grid nodes inside circle of interest (radiusMax)
xSourceInCirc = xSource(xSource.^2 + ySource.^2 <= radiusMax^2);
ySourceInCirc = ySource(xSource.^2 + ySource.^2 <= radiusMax^2);
zSourceInCirc = zSource(xSource.^2 + ySource.^2 <= radiusMax^2);

%Backp-project the rotated hologram to the surface of the source
SourceParameters.xGrid = xSourceInCirc;
SourceParameters.yGrid = ySourceInCirc;
SourceParameters.zGrid = zSourceInCirc;
SourceParameters.input = [];

FieldParameters.xGrid = xHoloRotated;
FieldParameters.yGrid = yHoloRotated;
FieldParameters.zGrid = zHoloRotated;
FieldParameters.input = holoRotated;
FieldParameters.dx = HologramSf.dx;
FieldParameters.dy = HologramSf.dy;

regime = 2;  % 2 Back-projection: P on a plane --> V on a sphere

[ vSourceInCirc ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

%Reshape the output velocity back at the rectangular grid for subsequent visualization
vSource = zeros(size(xSource));
vSource(xSource.^2 + ySource.^2 <= radiusMax^2) = vSourceInCirc;

%GUI part to change the radius of curvature, re-backporapagate the field if desired, and save the results
SourceParameters.xGrid = xSource;
SourceParameters.yGrid = ySource;
SourceParameters.zGrid = zSource;
SourceParameters.input = [];

if exist('xZeroPhase','var') ~= 1
  xZeroPhase = [];
end

if exist('yZeroPhase','var') ~= 1
  yZeroPhase = [];
end

gui_find_roc(vSource, expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature, xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, xZeroPhase, yZeroPhase, radiusMax, dxSource, dySource);

disp('xDDx simulation completed!');