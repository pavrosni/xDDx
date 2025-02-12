% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% ***Automatic Single-Frequency Hologram Rotation and Backprojection***
%
% This tool rotates and shifts a measured single-frequency hologram so that
% the aligned hologram is perpendicular to the axis of symmetry of the 
% transducer and back-project it to its surface. The script reads 
% predetermined position parameters 'alignmentParametersFileName' from a 
% file. This file can be either saved using the Automatic Alignment Tool 
% (bp_auto_alignment_sf.m) or constructed manually
% 
% INPUT DATA FORMAT FOR A SINGLE-FREQUENCY HOLOGRAM: a MAT or XLSX file 
% 'inputParametersFileName' with the following variables
% 'Geometry' struct with fields:
%     'radiusOfCurvature': nominal radius of curvature of the transducer in m
%     'apertureMin': characteristic minimum aperture of the transducer in m     
%     'apertureMax': characteristic maximum aperture of the transducer in m
%                    (apertureMax = apertureMin for circular transducers)     
%
% 'Medium' struct with fields:
%     'soundSpeed': sound speed in m/s
%     'density': density in kg/m^3
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
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

libraryDir = '..\..\..\xDDx_lib'; %lib directory

%I. SIMULATION REGIME
simulationDevice = 'cuda'; %simulation device: 'cuda' or 'cpu'


%II. INPUT PARAMETERS
inputParametersFileName = '..\..\data_for_examples\spherical\holo_data_spherical_sf.mat'; %input data with transducer, hologram, and medium parameters in '.mat' or '.xlsx' format
alignmentParametersFileName = '..\..\data_for_examples\spherical\alig_param_spherical.mat'; %data structure with position parameters (can be saved using 'bp_auto_alignment_sf.m')


%III. SOURCE PARAMETERS
nxSource = 201; %number of points for the Cartesian grid of the source along the x-dimension
nySource = 201; %number of points for the Cartesian grid of the source along the y-dimension
dxSource = 0.5e-3; %x grid step in m
dySource = 0.5e-3; %y grid step in m
xZeroPhase = -25e-3; %(optional parameter) x-coordinate in m of the zero-phase point on the surface of the source
yZeroPhase = -25e-3; %(optional parameter) y-coordinate in m of the zero-phase point on the surface of the source

%VI. TECHNICAL PARAMETERS
minHoloShiftInWl = 4; %minimum shift in wavelengths of the central frequency for the rotated hologram
apertureReserve = 0.1; %fraction of aperture used to increase the nominal output window
saveTransducer = false; %true to save the boundary condition for the transducer at HologramSf.frequency as a 'TransducerSf' structure

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; % cascade shift for each new figure in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(libraryDir));

%Load input data
[~, ~, inputExt] = fileparts(inputParametersFileName);
if strcmpi(inputExt, '.mat')
    load(inputParametersFileName);
else
    [Geometry, HologramSf, Medium] = read_hologram_sf_from_xls(inputParametersFileName);
end

load(alignmentParametersFileName);

%Calculate supplementary parameters
radiusMax = 0.5*(1+apertureReserve)*max(Geometry.apertureMin, Geometry.apertureMax);

%Calculate the holo-shift in m
frequency = HologramSf.frequency;
minHoloShift = minHoloShiftInWl*Medium.soundSpeed/frequency;

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

%Forward-project and rotate hologram
xHolo = HologramSf.xGrid;
yHolo = HologramSf.yGrid;
zHolo = zeros(size(xHolo));

expSign = HologramSf.expSign;

SourceParameters.xGrid = xHolo;
SourceParameters.yGrid = yHolo;
SourceParameters.zGrid = zHolo;
SourceParameters.dx = HologramSf.dx;
SourceParameters.dy = HologramSf.dy;
SourceParameters.input = HologramSf.complexPressureAmplitude;

regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points

isTransient = false; %set to true for the transient regime, false for the single-frequency regime

[xHoloRotated, yHoloRotated, zHoloRotated, holoRotated] = rotate_hologram(xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, Medium, ServiceParameters, zMaxAcoustCoord, minHoloShift);

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

%Back-project the rotated hologram to the surface of the source
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

if isSphericalSource
    regime = 2; % 2 Back-projection: P on a plane --> V on a sphere
    
    %Rayleigh simulator function with the complex velocity amplitude output
    [ vSourceInCirc ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);
else
    regime = 1; % 1 Back-projection: P on a plane --> V on a plane
    
    %Rayleigh simulator function with the complex velocity amplitude output
    [ vSourceInCirc ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);
end

%Reshape the output velocity back at the rectangular grid for subsequent visualization 
vSource = zeros(size(xSource));
vSource(xSource.^2 + ySource.^2 <= radiusMax^2) = vSourceInCirc;

%Normalize phase given the predetermined zero-phase point (xZeroPhase, yZeroPhase)
if (exist('xZeroPhase','var') ~= 0) && (exist('yZeroPhase','var') ~= 0)
    [~, iZeroPhase] = min((xSource(:) - xZeroPhase).^2 + (ySource(:) - yZeroPhase).^2);
    angSource = shift_2pi(angle(vSource) - angle(vSource(iZeroPhase)));
    angSource(xSource.^2 + ySource.^2 > radiusMax^2) = 4;
else
    angSource = angle(vSource);
end

%Plot complex velocity amplitude and phase at the surface of the transducer
positionCurrent = get(groot,'DefaultFigurePosition');

figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,abs(vSource)*1e3);
set(gca,'YDir','normal');
colormap('jet');
axis equal;
axis tight;
title('Vibrational velocity amplitude at the array surface, mm/s');
xlabel('{\itx}, mm');
ylabel('{\ity}, mm');
if isOctave
  xlim(xSource([1 end])*1e3)
  ylim(ySource([1 end])*1e3)
end
colorbar;

figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,angSource);
set(gca,'YDir','normal');
colormap('hsv');
axis equal;
axis tight;
title('Vibrational velocity phase at the array surface, rad');
xlabel('{\itx}, mm');
ylabel('{\ity}, mm');
if isOctave
  xlim(xSource([1 end])*1e3)
  ylim(ySource([1 end])*1e3)
end
colorbar;

if saveTransducer

    if ~isSphericalSource
        radiusOfCurvature = [];
    end
    save_transducer(isSphericalSource, expSign, HologramSf.frequency, xSource, ySource, zSource, radiusOfCurvature, dxSource, dySource, vSource);

end

disp('xDDx simulation completed!');