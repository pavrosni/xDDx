% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% ***Simple Hologram Forwardpropagation***
%
% Forward-projects a measured single-frequency hologram to obtain a 3D, 2D,
% 1D, or 0D field. The output acoustic pressure amplitude is plotted in 3D,
% 2D, 1D, or 0D with an option of a slice-by-slice representation in 3D case
%
%
% INPUT DATA FORMAT FOR A SINGLE-FREQUENCY HOLOGRAM: a MAT or XLSX file
% 'inputParametersFileName' with the following variables 
%
% 'Geometry' struct with fields
%     'radiusOfCurvature': nominal radius of curvature of the transducer in m 
%     'apertureMin': characteristic minimum aperture of the transducer in m     
%     'apertureMax': characteristic maximum aperture of the transducer in m
%                    (apertureMax = apertureMin for circular transducers)     
%
% 'Medium' struct with fields
%     'soundSpeed': sound speed in m/s
%     'density': density in kg/m^3
%
% 'HologramSf' struct with fields
%     'expSign': +1 or -1, depending on the exponent sign convention exp(+
%     1i * omega * t) or exp(- 1i * omega * t). 
% 	             E.g., the "fft" function in MATLAB utilizes the  exp(+ 1i *
% 	             omega * t) convention, so in this case, expSign is +1 
%     'frequency': the operating frequency of the transducer in Hz
%     'xGrid': vector or matrix with the x-coordinates in m at each grid
%     node of the hologram  (Cartesian grid only)
%     'yGrid': vector or matrix with the y-coordinates in m at
%     corresponding grid nodes of the hologram  (Cartesian grid only)
%     'dx': x-step of the hologram Cartesian grid in m
%     'dy': y-step of the hologram Cartesian grid in m
%     'complexPressureAmplitude': vector or matrix with the complex
%     pressure amplitude values in Pa at corresponding grid nodes of the hologram 
%     'zPosition': approximate z-position in m of the hologram, assuming
%     zero at the apex of the transducer 
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

%I. SIMULATION DEVICE
simulationDevice = 'cuda'; %simulation device: 'cuda' or 'cpu'

%II. INPUT PARAMETERS
inputParametersFileName = '..\..\data_for_examples\spherical\holo_data_spherical_sf.mat'; %input data with transducer, single-frequency hologram, and medium parameters in '.mat' or '.xlsx' format

%III. RECTANGULAR WINDOW FOR FIELD SIMULATION
xFieldBegin = -10e-3; %x, y, and z limits in m for the rectangular simulation region. Set "Begin" equal to "End" for specific coordinates to reduce the number of dimensions.
xFieldEnd   =  10e-3;
yFieldBegin = -10e-3;
yFieldEnd   = 10e-3;
zFieldBegin = 70e-3;
zFieldEnd   = 110e-3;

dxField = 0.25e-3; %x, y, and z grid step in m for the rectangular simulation region. Steps for singleton dimensions (where "Begin" = "End") are ignored and can be omitted.
dyField = 0.25e-3;
dzField = 0.5e-3;

%IV. ISOLEVEL PARAMETERS FOR 3D (ignored or can be omitted in 2D/1D/0D cases)
levelArray = [0.8 0.5 0.1]; %isolevels related to the pressure maximum for extracting the isosurface of the simulated 3D field
transparencyArray = [0.6 0.5 0.1]; %transparency of the isosurfaces for the levels from levelArray


%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
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

%Set geometric parameters
xHolo = HologramSf.xGrid;
yHolo = HologramSf.yGrid;
zHolo = HologramSf.zPosition * ones(size(xHolo));

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

%Forwardpropagate the hologram to calculate the complex pressure amplitude at the nodes of the grid
SourceParameters = [];
SourceParameters.xGrid = xHolo;
SourceParameters.yGrid = yHolo;
SourceParameters.zGrid = zHolo;
SourceParameters.dx = HologramSf.dx;
SourceParameters.dy = HologramSf.dy;
SourceParameters.input = HologramSf.complexPressureAmplitude;

FieldParameters = [];
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;

regime = 5; % 5 Forwardpropagation: P on a plane --> P at an arbitrary set of points

isTransient = false; %set to true for the transient regime, false for the single-frequency regime

%Rayleigh simulator function with the complex pressure amplitude output
[ pField3D ] = rayleigh_simulator(HologramSf.expSign, HologramSf.frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

%GUI part to show 3D/2D/1D/0D results in a slice-by-slice format and save them

if abs(ndims(squeeze(pField3D)) - 3) < eps
    gui_plot_3d(xField3D, yField3D, zField3D, pField3D, levelArray, transparencyArray);
elseif ~isvector(squeeze(pField3D))
    [xField2D, yField2D, zField2D, pField2D] = plot_2d_field(xField3D, yField3D, zField3D, pField3D);
elseif numel(squeeze(pField3D)) > 1
    [xField1D, yField1D, zField1D, pField1D] = plot_1d_field(xField3D, yField3D, zField3D, pField3D);
else
    [xField0D, yField0D, zField0D, pField0D] = disp_0d_field(xField3D, yField3D, zField3D, pField3D);
end

disp('xDDx simulation completed!');