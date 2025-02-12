% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% An example of a simulation script for the numerical simulation of the 3D,
% 2D, 1D, or 0D field field of a user-defined transducer operating in a
% single-frequency regime. The output acoustic pressure amplitude is plotted in
% 3D,2D, 1D, or 0D with an option of a slice-by-slice representation in 3D
% case. Transducer parameters are given by the 'TransducerSf' structure.
%
% This example uses a multi-element transducer with a nonregular
% vibrational velocity pattern: a 109-element 1-MHz spherically shaped
% array with a rectangular aperture, an oval opening and a randomized fully
% populated pattern of elements (see design method in
% [Rosnitskiy et al., DOI 10.1109/TUFFC.2018.2800160]). Phases at the array
% elements were specifically set to provide a 10-mm electronic steering of
% the focus in the x-direction. Three different examples of transducer's
% grid steps are presented.
%
% Note that you can use structures 'TransducerSf' saved from the GUIs in
% corresponding tools:
%
% 'xDDx\examples\alignment_tools\bp_auto_alignment_sf.m'
% 'xDDx\examples\alignment_tools\bp_predefined_alignment_transient.m'
% (button 'Save TransducerSf')
%
% The output of the script shows the vibrational velocity amplitude/phase
% at the surface of the transducer, the 3D field and the transducer surface,
% and a GUI window with 3D field simulation results in a slice-by-slice
% format.
%
% INPUT DATA FORMAT WITH TRANSDUCER PARAMETERS OPERATING IN A
% SINGLE-FREQUENCY REGIME: a MAT or XLSX file
% 'inputParametersFileName' with the following variable
%
% 'TransducerSf' struct with fields:
%   'expSign': +1 or -1, depending on the exponent sign convention
%       exp(+ 1i * omega * t) or exp(- 1i * omega * t). E.g., the "fft"
%       function in MATLAB utilizes the  exp(+ 1i * omega * t) convention,
%       so in this case, expSign is +1
%   'frequency': the operating frequency of the transducer in Hz
%   'radiusOfCurvature': (optional) radius of curvature of the transducer in m
%   'xGrid': vector or matrix with the x-coordinates in m at each grid node
%       of the transducer surface (Cartesian grid only)
%   'yGrid': vector or matrix with the y-coordinates in m at each grid node
%       of the transducer surface (Cartesian grid only)
%   'zGrid': vector or matrix with the z-coordinates in m at each grid node
%       of the transducer surface (Cartesian grid only). Users can omit this
%       variable, in such cases, the toolbox will automatically determine the
%       transducer type - spherical or flat.
%       For a spherical transducer, zGrid nodes will be set according to
%       the sphere equation:
%       zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xGrid.^2 - yGrid.^2)
%       For a flat transducer, all zGrid nodes will be set to zero:
%       zGrid = zeros(size(xGrid)).
%   'dx': x-step of the transducer Cartesian grid in m
%   'dy': y-step of the transducer Cartesian grid in m
%   'complexVelocityAmplitude': complex velocity amplitude in m/s at each
%       node of the transducer surface grid
%
%  XLSX template path is
%  'xDDx\examples\data_for_examples\xlsx_templates\transducer_sf.xls'
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
% inputParametersFileName =
% '..\data_for_examples\simulation\fpa_steering_10mm_fine_grid.mat'; %input with the transducer data in 'mat' or 'xlsx' format
inputParametersFileName = '..\data_for_examples\simulation\fpa_steering_10mm_standard_grid.mat'; %input with the transducer data in '.mat' or '.xlsx' format
% inputParametersFileName = '..\data_for_examples\simulation\fpa_steering_10mm_coarse_grid.mat'; %input with the transducer data in '.mat' or .'xlsx' format
% inputParametersFileName = '..\data_for_examples\simulation\fpa_steering_10mm_standard_grid.xlsx'; %input with the transducer data in '.mat' or '.xlsx' format

%IV. 3D WINDOW FOR FIELD SIMULATION
xFieldBegin = -15e-3; %x, y, and z limits in m for the rectangular simulation region. Set "Begin" equal to "End" for specific coordinates to reduce the number of dimensions.
xFieldEnd   =  15e-3;
yFieldBegin = -5e-3;
yFieldEnd   =  5e-3;
zFieldBegin = 60e-3;
zFieldEnd   = 100e-3;

dxField = 0.25e-3;  %x, y, and z grid step in m for the rectangular simulation region. Steps for singleton dimensions (where "Begin" = "End") are ignored and can be omitted.
dyField = 0.25e-3;
dzField = 0.5e-3;


%V. ISOLEVEL PARAMETERS
levelArray = [0.8 0.5 0.3]; %isolevels related to the pressure maximum for extracting the isosurface of the simulated 3D field
transparencyArray = [0.6 0.5 0.2]; %transparency of the isosurfaces for the levels from levelArray


%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; %cascade shift for each new figure in pixels
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
    TransducerSf = read_transducer_sf_from_xls(inputParametersFileName);
end

[TransducerSf] = check_transducer_z(TransducerSf);

if isvector(TransducerSf.xGrid)
    [TransducerSf] = reshape_transducer_dim(TransducerSf);
end

isTxMeshgrid = is_meshgrid(TransducerSf.xGrid, TransducerSf.yGrid);
meshgridWarning = 'Transducer surface is not shown because your input grid (xGrid, yGrid) is not in meshgrid format!';

%Set geometric parameters

expSign = TransducerSf.expSign;
frequency = TransducerSf.frequency;

if isfield(TransducerSf, 'radiusOfCurvature')
    if ~isempty(TransducerSf.radiusOfCurvature)
        radiusOfCurvature = TransducerSf.radiusOfCurvature;
    end
end
isSphericalSource = exist('radiusOfCurvature', 'var') ~= 0;

xSource = TransducerSf.xGrid;
ySource = TransducerSf.yGrid;
zSource = TransducerSf.zGrid;

xFieldVector = xFieldBegin;
yFieldVector = yFieldBegin;
zFieldVector = zFieldBegin;

if abs(xFieldEnd - xFieldBegin) > eps
    xFieldVector = xFieldBegin : dxField : xFieldEnd;
end

if abs(yFieldEnd - yFieldBegin) > eps
    yFieldVector = yFieldBegin : dyField : yFieldEnd;
end

if abs(zFieldEnd - zFieldBegin) > eps
    zFieldVector = zFieldBegin : dzField : zFieldEnd;
end

%Generate 3D grid for simulation output
[xField3D, yField3D, zField3D] = meshgrid(xFieldVector, yFieldVector, zFieldVector);

%Forward-project the hologram to calculate the complex pressure amplitude at the nodes of the grid

withinTheActiveSurface = (abs(TransducerSf.complexVelocityAmplitude)>eps); %surface mask of non-zero velocity

SourceParameters = [];
SourceParameters.xGrid = xSource(withinTheActiveSurface);
SourceParameters.yGrid = ySource(withinTheActiveSurface);
SourceParameters.zGrid = zSource(withinTheActiveSurface);
SourceParameters.dx = TransducerSf.dx;
SourceParameters.dy = TransducerSf.dy;
SourceParameters.input = TransducerSf.complexVelocityAmplitude(withinTheActiveSurface);

FieldParameters = [];
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;

isTransient = false; %set to true for the transient regime, false for the single-frequency regime

if isSphericalSource
    regime = 3; % 3 Forward-projection: V on a sphere --> P at an arbitrary set of points
    %Rayleigh simulator function with the complex pressure amplitude output
    [ pField3D ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

else

    regime = 4; % 4 Forward-projection: V on a plane --> P at an arbitrary set of points
    %Rayleigh simulator function with the complex pressure amplitude output
    [ pField3D ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

end

%Find the pressure maximum point (xMax3D, yMax3D, zMax3D)
[pMax3D,iMax3D] = max(abs(pField3D(:)));
xMax3D = xField3D(iMax3D);
yMax3D = yField3D(iMax3D);
zMax3D = zField3D(iMax3D);

positionCurrent = get(groot,'DefaultFigurePosition');
if isTxMeshgrid

    %Plot complex velocity amplitude and phase at the surface of the transducer
    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    hvAmpl = imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,abs(TransducerSf.complexVelocityAmplitude)*1e3);
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
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,angle(TransducerSf.complexVelocityAmplitude));
    set(gca,'YDir','normal');
    colormap('hsv');
    caxis([-pi pi]);
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

else

    warning(meshgridWarning);

end


if abs(ndims(squeeze(pField3D)) - 3) < eps
    %Plot 3D field and source surface
    figure;
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

    if isTxMeshgrid

        cDatavAmpl = get(hvAmpl, 'CData');

        %Normalize the CData to the range of the colormap
        cDataNormalized = (cDatavAmpl - min(cDatavAmpl(:))) / (max(cDatavAmpl(:)) - min(cDatavAmpl(:)));
        cDataNormalized = round(cDataNormalized * (size(cmap, 1) - 1)) + 1;

        %Convert to RGB using the colormap
        rgbImageAmpl = ind2rgb(cDataNormalized, cmap);

        activeSurfaceMask = ones(size(TransducerSf.complexVelocityAmplitude));
        activeSurfaceMask(abs(TransducerSf.complexVelocityAmplitude)<eps)=0;

        if ~isOctave
            h = surface(TransducerSf.xGrid*1e3, TransducerSf.yGrid*1e3, TransducerSf.zGrid*1e3, rgbImageAmpl, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
            set(h, 'AlphaData', activeSurfaceMask, 'FaceAlpha', 'texturemap', 'CDataMapping', 'direct');
            camlight;
        else
            surface(TransducerSf.xGrid*1e3, TransducerSf.yGrid*1e3, TransducerSf.zGrid*1e3, rgbImageAmpl, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
            zlim([0 max(zField3D(:))*1e3]);
        end
        axis equal;
        if ~isOctave
            set(h, 'AlphaData', activeSurfaceMask, 'FaceAlpha', 'texturemap', 'CDataMapping', 'direct');
            view(3);
            camlight;
        else
            zlim([0 max(zField3D(:))*1e3]);
            view(3);
        end

    else

        if ~isOctave
            view(3);
            camlight;
        else
            zlim([0 max(zField3D(:))*1e3]);
            view(3);
        end

    end

    gui_plot_3d(xField3D,yField3D, zField3D, pField3D, levelArray, transparencyArray, false);
elseif ~isvector(squeeze(pField3D))
    [xField2D, yField2D, zField2D, pField2D] = plot_2d_field(xField3D, yField3D, zField3D, pField3D);
elseif numel(squeeze(pField3D)) > 1
    [xField1D, yField1D, zField1D, pField1D] = plot_1d_field(xField3D, yField3D, zField3D, pField3D);
else
    [xField0D, yField0D, zField0D, pField0D] = disp_0d_field(xField3D, yField3D, zField3D, pField3D);
end


disp('xDDx simulation completed!');
