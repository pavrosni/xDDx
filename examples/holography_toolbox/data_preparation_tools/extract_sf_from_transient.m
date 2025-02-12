% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%***Extract a single-frequency hologram from a transient hologram***
% 
% Reads transient hologram data from a file, specified by 
% inputParametersFileName, and extracts a single-frequency component at a 
% frequency of frequencyNominal. The output frequency is as close to the 
% desired frequencyNominal as the spectral sampling of the discrete signal 
% allows. The result is saved as outFilename.
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
% OUTPUT DATA FORMAT FOR A SINGLE-FREQUENCY HOLOGRAM: a MAT file
% 'outFilename' with the following variables
%
% 'Geometry' struct with fields:
%   'radiusOfCurvature': (optional) nominal radius of curvature of the 
%                        transducer in m. Omit this variable or define it 
%                        as an empty value [] if your transducer is flat.
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

%II. OUTPUT FILE NAME
outFilename = 'holo_data_sf.mat'; %output data with transducer, single-frequency hologram, and medium parameters

%III. TRANSDUCER PARAMETERS
frequencyNominal = 1.25e6; %desired frequency of the output single-frequency hologram in Hz 
apertureMin = 87e-3; % characteristic minimum aperture of the transducer in m
apertureMax = 87e-3; % characteristic maximum aperture of the transducer in m (apertureMax = apertureMin for circular transducers)
radiusOfCurvature = 87e-3; %(optional parameter) nominal radius of curvature of the transducer in m. Omit this variable if your transducer is flat. 

%TECHNICAL PARAMETERS
expSign = +1; %+1 or -1, depending on the exponent sign convention exp(+ 1i * omega * t) or exp(- 1i * omega * t), +1 by default. 

%SERVICE PARAMETERS
figuresCascadeShift = 50; %cascade shift for each new figure in pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(libraryDir));

%Load input data
load(inputParametersFileName);

needToReshape = isvector(HologramTr.xGrid);
if needToReshape
    [HologramTr] = reshape_transducer_dim(HologramTr);
end

[HologramSf, frequency] = extract_sf(HologramTr,expSign,frequencyNominal);


%Create Geometry structure
Geometry = [];

if exist('radiusOfCurvature','var') ~= 0
    Geometry.radiusOfCurvature = radiusOfCurvature;
end

Geometry.apertureMin = apertureMin;
Geometry.apertureMax = apertureMax;

if needToReshape
    HologramSf.xGrid = squeeze(HologramSf.xGrid);
    HologramSf.yGrid = squeeze(HologramSf.yGrid);
    HologramSf.complexPressureAmplitude = squeeze(HologramSf.complexPressureAmplitude);
end

%Save output data file
save(outFilename, 'Geometry', 'HologramSf', 'Medium');

%Show results
positionCurrent = get(groot,'DefaultFigurePosition');

if ~needToReshape
    %Plot amplitude
    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    imagesc(abs(HologramSf.complexPressureAmplitude)*1e-3);
    colormap jet;
    colorbar;
    xlabel('# of x-sample');
    ylabel('# of y-sample');
    title(['Pressure amplitude hologram [kPa] at ' num2str(frequency*1e-6,5) ' MHz' ]);
    axis equal;
    axis tight;

    %Plot phase
    figure;
    set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
    positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
    imagesc(angle(HologramSf.complexPressureAmplitude));
    colormap hsv;
    colorbar;
    xlabel('# of x-sample');
    ylabel('# of y-sample');
    title(['Pressure phase hologram [rad.] at ' num2str(frequency*1e-6,5) ' MHz' ]);
    axis equal;
    axis tight;

else

  disp('HologramSf is saved!');  

end

disp('xDDx simulation completed!');