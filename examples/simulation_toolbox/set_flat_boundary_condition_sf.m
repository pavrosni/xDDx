% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% This tool may be useful for back-projecting a boundary condition from a 
% spherical surface to a flat surface can be used to set a boundary 
% condition for spherically shaped transducers on a plane. It can serve as 
% one of the possible simple methods for setting flat boundary conditions 
% for spherically shaped transducers in other simulation software such as 
% k-Wave, FOCUS, and mSOUND.
%
% Transducer parameters are provided by a ‘TransducerSf’ structure. 
% The script back-projects the vibrational velocity at the surface of the 
% transducer to a plane at the apex of the transducer and then numerically 
% simulates the 3D field based on the obtained flat boundary condition to 
% test its accuracy. 
%
% This example uses a multi-element transducer with a nonregular
% vibrational velocity pattern: a 109-element 1-MHz spherically shaped
% array with a rectangular aperture, an oval opening and a randomized fully
% populated pattern of elements (see design method in
% [Rosnitskiy et al., DOI 10.1109/TUFFC.2018.2800160]). Phases at the array
% elements were specifically set to provide a 10-mm electronic steering of
% the focus in the x-direction. Three different examples of transducer’s
% grid steps are presented.
%
% The output of the script shows the vibrational velocity amplitude/phase
% at the boundary condition plane, the 3D field,  and a GUI window with 3D
% field simulation results in a slice-by-slice format.
%
% Note that you can use structures 'TransducerSf' saved from the GUIs in
% corresponding tools:
%
% 'xDDx\examples\alignment_tools\bp_auto_alignment_sf.m' 
% 'xDDx\examples\alignment_tools\bp_predefined_alignment_transient.m' 
% (button 'Save TransducerSf') 
%
% INPUT DATA FORMAT WITH TRANSDUCER PARAMETERS OPERATING IN A
% SINGLE-FREQUENCY REGIME: a MAT file
% 'inputParametersFileName' with the following variable
%
% 'TransducerSf' struct with fields:
%   'expSign': +1 or -1, depending on the exponent sign convention
%       exp(+ 1i * omega * t) or exp(- 1i * omega * t). E.g., the "fft"
%       function in MATLAB utilizes the  exp(+ 1i * omega * t) convention,
%       so in this case, expSign is +1
%   'frequency': the operating frequency of the transducer in Hz
%   'radiusOfCurvature': radius of curvature of the transducer in m
%   'xGrid': vector or matrix with the x-coordinates in m at each grid node
%       of the transducer surface (Cartesian grid only)
%   'yGrid': vector or matrix with the y-coordinates in m at each grid node
%       of the transducer surface (Cartesian grid only)
%   'zGrid': vector or matrix with the z-coordinates in m at each grid node
%       of the transducer surface (Cartesian grid only)
%   'dx': x-step of the transducer Cartesian grid in m
%   'dy': y-step of the transducer Cartesian grid in m
%   'complexVelocityAmplitude': complex velocity amplitude in m/s at each
%       node of the transducer surface grid
%
% Here, x and y are the transverse coordinates of the transducer, and z is 
% the coordinate along the direction of the beam. See details in the manual.

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%EDITABLE CODE%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

libraryDir = '..\..\xDDx_lib'; %lib directory

%I. SIMULATION DEVICE
simulationDevice = 'cuda'; %simulation device: 'cuda' or 'cpu'


%II. MEDIUM PARAMETERS
Medium.soundSpeed = 1500; %sound speed in m/s
Medium.density = 1000; %density in kg/m^3


%III. INPUT PARAMETERS
% inputParametersFileName = '..\data_for_examples\simulation\fpa_steering_10mm_fine_grid.mat'; %input with the transducer data
inputParametersFileName = '..\data_for_examples\simulation\fpa_steering_10mm_standard_grid.mat'; %input with the transducer data
% inputParametersFileName = '..\data_for_examples\simulation\fpa_steering_10mm_coarse_grid.mat'; %input with the transducer data


%FLAT BOUNDARY CONDITION
nxFlatBoundary = 181; %number of points for the Cartesian grid of the flat boundary condition along the x-dimension
nyFlatBoundary = 141; %number of points for the Cartesian grid of the flat boundary condition along the y-dimension
dxFlatBoundary = 0.7e-3; %x and y grid step in m for the flat boundary condition
dyFlatBoundary = 0.7e-3;


%IV. 3D WINDOW FOR FIELD SIMULATION
xFieldBegin = -15e-3; %x, y, and z limits in m for the rectangular 3D simulation region
xFieldEnd   =  15e-3;
yFieldBegin = -5e-3;
yFieldEnd   =  5e-3;
zFieldBegin = 60e-3;
zFieldEnd   = 100e-3;

dxField = 0.25e-3; %x, y, and z grid step in m for the rectangular 3D simulation region
dyField = 0.25e-3;
dzField = 0.5e-3;


%V. ISOLEVEL PARAMETERS
levelArray = [0.8 0.5 0.3]; %isolevels related to the pressure maximum for extracting the isosurface of the simulated 3D field
transparencyArray = [0.6 0.5 0.2]; %transparency of the isosurfaces for the levels from levelArray

%TECHNICAL PARAMETERS
shiftHoloDistInWl = 10; %minimum back-projection distance in wavelengths for the flat boundary condition

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; %cascade shift for each new figure in pixels
saveBoundaryCondition = false; %true to save the flat boundary condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

addpath(genpath(libraryDir));

%Load input data
load(inputParametersFileName);

if isvector(TransducerSf.xGrid)
    [TransducerSf] = reshape_transducer_dim(TransducerSf);
end

%Set geometric parameters
expSign = TransducerSf.expSign;
frequency = TransducerSf.frequency;
if isfield(TransducerSf, 'radiusOfCurvature')
    if ~isempty(TransducerSf.radiusOfCurvature)
        radiusOfCurvature = TransducerSf.radiusOfCurvature;
    end
end
isSphericalSource = exist('radiusOfCurvature', 'var') ~= 0;
if isSphericalSource
  isSphericalSource = isSphericalSource*(~isempty(radiusOfCurvature));
end

if ~isSphericalSource
  error('Your transducer is already flat! There is no need to apply the boundary condition transfer.');
end

%Calculate the shift distance in m
shiftHoloDist = shiftHoloDistInWl*Medium.soundSpeed/frequency;

%Build the boundary condition grid that includes the center of symmetry of the transducer
[xFlatBoundary, yFlatBoundary] = build_flat_grid_centered(nxFlatBoundary,dxFlatBoundary,nyFlatBoundary,dyFlatBoundary);
zFlatBoundary = zeros(size(xFlatBoundary));

%Back-project the source vibrational velocity to calculate the complex
% pressure amplitude at the nodes of the flat boundary condition at z = 0.
% If the source is too close to the point z = 0, additional back projection
% and forward projection are done to fulfill the predefined minimum 
% back-projection distance shiftHoloDistInWl
if abs(min(TransducerSf.zGrid(:))) < shiftHoloDist
    
    %Save parameters of the desired boundary condition at z = 0
    xFlatBoundaryDesired = xFlatBoundary;
    yFlatBoundaryDesired = yFlatBoundary;
    zFlatBoundaryDesired = zFlatBoundary;
    
    %Increase the size of the shifted flat boundary condition as it is more
    %distant from the focus than the desired one
    maxTransverceSize = max(nxFlatBoundary*dxFlatBoundary,nyFlatBoundary*dyFlatBoundary);
    rezervedSize = 2*shiftHoloDist*(0.5*maxTransverceSize/TransducerSf.radiusOfCurvature);
    nxFlatBoundary = nxFlatBoundary + fix(rezervedSize/dxFlatBoundary) + 1;
    nyFlatBoundary = nyFlatBoundary + fix(rezervedSize/dyFlatBoundary) + 1;
    [xFlatBoundary, yFlatBoundary] = build_flat_grid_centered(nxFlatBoundary,dxFlatBoundary,nyFlatBoundary,dyFlatBoundary);
    zFlatBoundary = -shiftHoloDist*ones(size(xFlatBoundary));

    %Back-project the source vibrational velocity to calculate the complex
    %pressure amplitude at the nodes of the shifted flat boundary condition
    SourceParameters = [];
    SourceParameters.xGrid = xFlatBoundary;
    SourceParameters.yGrid = yFlatBoundary;
    SourceParameters.zGrid = zFlatBoundary;
    FieldParameters = [];
    withinTheActiveSurface = (abs(TransducerSf.complexVelocityAmplitude)>eps);
    FieldParameters.xGrid = TransducerSf.xGrid(withinTheActiveSurface);
    FieldParameters.yGrid = TransducerSf.yGrid(withinTheActiveSurface);
    FieldParameters.zGrid = TransducerSf.zGrid(withinTheActiveSurface);
    FieldParameters.dx = TransducerSf.dx;
    FieldParameters.dy = TransducerSf.dy;
    FieldParameters.input = TransducerSf.complexVelocityAmplitude(withinTheActiveSurface);
    regime = 6; % 6 Back-projection: V on a sphere --> P on a planes
    isTransient = false;
    [ pFlatBoundary] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);


    %Forward-project the shifted flat boundary condition to calculate the
    %complex pressure amplitude at the nodes of the desired flat boundary
    %condition
    SourceParameters = [];
    SourceParameters.xGrid = xFlatBoundary;
    SourceParameters.yGrid = yFlatBoundary;
    SourceParameters.zGrid = zFlatBoundary;
    SourceParameters.dx = dxFlatBoundary;
    SourceParameters.dy = dyFlatBoundary;
    SourceParameters.input = pFlatBoundary;
    FieldParameters = [];
    FieldParameters.xGrid = xFlatBoundaryDesired;
    FieldParameters.yGrid = yFlatBoundaryDesired;
    FieldParameters.zGrid = zFlatBoundaryDesired;
    regime = 5; % 5 Forward-projection: P on a plane --> P at an arbitrary set of points
    isTransient = false;
    [ pFlatBoundaryDesired] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

    %Go back to the desired boundary condition at z = 0
    pFlatBoundary = pFlatBoundaryDesired;
    xFlatBoundary = xFlatBoundaryDesired;
    yFlatBoundary = yFlatBoundaryDesired;
    zFlatBoundary = zFlatBoundaryDesired;

else

    %Back-project the source vibrational velocity to calculate the complex
    %pressure amplitude at the nodes of the flat boundary condition
    SourceParameters = [];
    SourceParameters.xGrid = xFlatBoundary;
    SourceParameters.yGrid = yFlatBoundary;
    SourceParameters.zGrid = zFlatBoundary;
    FieldParameters = [];
    withinTheActiveSurface = (abs(TransducerSf.complexVelocityAmplitude)>eps);
    FieldParameters.xGrid = TransducerSf.xGrid(withinTheActiveSurface);
    FieldParameters.yGrid = TransducerSf.yGrid(withinTheActiveSurface);
    FieldParameters.zGrid = TransducerSf.zGrid(withinTheActiveSurface);
    FieldParameters.dx = TransducerSf.dx;
    FieldParameters.dy = TransducerSf.dy;
    FieldParameters.input = TransducerSf.complexVelocityAmplitude(withinTheActiveSurface);
    regime = 6; % 6 Back-projection: V on a sphere --> P on a planes
    isTransient = false;
    [ pFlatBoundary] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);
end

%Use the angular spectrum method to convert the pressure boundary condition into the velocity one
vFlatBoundary = pressure_to_velocity_flat_surface(pFlatBoundary, dyFlatBoundary, dxFlatBoundary, frequency, Medium);


%Generate 3D grid for simulation output
xFieldVector = xFieldBegin : dxField : xFieldEnd;
yFieldVector = yFieldBegin : dyField : yFieldEnd;
zFieldVector = zFieldBegin : dzField : zFieldEnd;
[xField3D, yField3D, zField3D] = meshgrid(xFieldVector, yFieldVector, zFieldVector);

%Forward-project the flat boundary condition to calculate the complex pressure amplitude at the nodes of the 3D field grid
SourceParameters = [];
SourceParameters.xGrid = xFlatBoundary;
SourceParameters.yGrid = yFlatBoundary;
SourceParameters.zGrid = zFlatBoundary;
SourceParameters.dx = dxFlatBoundary;
SourceParameters.dy = dyFlatBoundary;
SourceParameters.input = vFlatBoundary;

FieldParameters = [];
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;

regime = 4; % 4 Forward-projection: V on a plane --> P at an arbitrary set of points

isTransient = false; %set to true for the transient regime, false for the single-frequency regime

%Rayleigh simulator function with the complex pressure amplitude output
[ pField3D ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

%Find the pressure maximum point (xMax3D, yMax3D, zMax3D)
[pMax3D,iMax3D] = max(abs(pField3D(:)));
xMax3D = xField3D(iMax3D);
yMax3D = yField3D(iMax3D);
zMax3D = zField3D(iMax3D);

%Plot complex velocity amplitude and phase at the surface of the transducer
positionCurrent = get(groot,'DefaultFigurePosition');
figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
hvAmpl = imagesc(xFlatBoundary([1 end])*1e3,yFlatBoundary([1 end])*1e3,abs(vFlatBoundary)*1e3);
set(gca,'YDir','normal');
cmvAmpl = colormap('jet');
axis equal;
axis tight;
title('Vibrational velocity amplitude at the flat boundary condition, mm/s');
xlabel('{\itx}, mm');
ylabel('{\ity}, mm');
if isOctave
    xlim(xFlatBoundary([1 end])*1e3)
    ylim(yFlatBoundary([1 end])*1e3)
end
colorbar;

figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
imagesc(xFlatBoundary([1 end])*1e3,yFlatBoundary([1 end])*1e3,angle(vFlatBoundary));
set(gca,'YDir','normal');
colormap('hsv');
axis equal;
axis tight;
title('Vibrational velocity phase at the flat boundary condition, rad');
xlabel('{\itx}, mm');
ylabel('{\ity}, mm');
if isOctave
    xlim(xFlatBoundary([1 end])*1e3)
    ylim(yFlatBoundary([1 end])*1e3)
end
colorbar;


%Plot 3D field and the boundary conditon
figure;
hold on;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
cmap = colormap('jet');
climits = [0 1];
scalar=levelArray;
scalarclamped = scalar;
scalarclamped(scalar < climits(1)) = climits(1);
scalarclamped(scalar > climits(2)) = climits(2);
colorsLevels = interp1(linspace(climits(1), climits(2), size(cmap, 1)), ...
    cmap, ...
    scalarclamped);
for i = 1:length(levelArray)
    p=patch(isosurface(xField3D*1e3, yField3D*1e3, zField3D*1e3, abs(pField3D)/max(abs(pField3D(:))), levelArray(i)));
    set(p,'FaceColor',squeeze(colorsLevels(i,:)),'EdgeColor','none','FaceAlpha',transparencyArray(i));
end
axis equal;
axis tight;
xlabel('x, mm');
ylabel('y, mm');
zlabel('z, mm');
colormap jet;
titleFor3D = ['Pressure amplitude max is ' num2str(pMax3D) ' Pa at ' '(' num2str(xMax3D*1e3) ', ' num2str(yMax3D*1e3) ', ' num2str(zMax3D*1e3) ') mm' ];
if ~isOctave
    cbh = colorbar;
    caxis([0 1]);
    set(cbh,'XTick',flip(levelArray));
else
    titleFor3D = {titleFor3D, ['pressure amplitude iso-levels: ' num2str(levelArray)]};
end
title(titleFor3D);
view(3);

cDatavAmpl = get(hvAmpl, 'CData');
% Normalize the CData to the range of the colormap
cDataNormalized = (cDatavAmpl - min(cDatavAmpl(:))) / (max(cDatavAmpl(:)) - min(cDatavAmpl(:)));
cDataNormalized = round(cDataNormalized * (size(cmvAmpl, 1) - 1)) + 1;
% Convert to RGB using the colormap
rgbImageAmpl = ind2rgb(cDataNormalized, cmvAmpl);
image(xFlatBoundary([1 end])*1e3,yFlatBoundary([1 end])*1e3,rgbImageAmpl);
set(gca,'YDir','normal');
axis equal;
if ~isOctave
    camlight;
else
    zlim([0 max(zField3D(:))*1e3]);
end

%GUI part to show 3D results in a slice-by-slice format and save them
gui_plot_3d(xField3D,yField3D, zField3D, pField3D, levelArray, transparencyArray, false);

%Save boundary condition 
if saveBoundaryCondition
    xFlatBoundary = xFlatBoundary;
    yFlatBoundary = yFlatBoundary;
    zFlatBoundary = zFlatBoundary;

    outFileName = ['flat_boundary_condition_' char(datetime('now','Format','d_MMM_y_HH_mm_ss'))];

    save(outFileName, 'xFlatBoundary', 'yFlatBoundary', 'zFlatBoundary', 'vFlatBoundary', 'TransducerSf', 'Medium');
end

disp('xDDx simulation completed!');