% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
%
% This is an example simulation script for the numerical simulation of the 
% field of a transducer operating in a transient regime. 
% 
% Transducer parameters are provided by the 'TransducerTr' structure, which 
% can be read from a specified file. The output includes frame-by-frame 
% GUI displaying simulated pressure patterns in 2D simulation regions at 
% different moments in time. An option to render results in video format 
% is available. For 0D, 1D, and 3D simulations, raw results can be saved by
% the user and visualized as needed.
%
% Note: transient simulations may require more processing time compared to
% single-frequency simulations.
%
% INPUT DATA FORMAT WITH TRANSDUCER PARAMETERS OPERATING IN A
% TRANSIENT REGIME: a MAT file
% 'inputParametersFileName' with the following variable
%
% 'TransducerTr' struct with fields:
%   'radiusOfCurvature': (optional, can be omitted) radius of curvature of 
%                        the transducer in m. Omit this field or set it as 
%                        an empty vector if your transducer is flat
%     'xGrid': vector or matrix with x-coordinates of the Cartesian grid in
%              m at each grid node of the transducer surface
%     'yGrid': vector or matrix with y-coordinates of the Cartesian grid in
%              m at each grid node of the transducer surface
%     'zGrid': vector or matrix with z-coordinates of the Cartesian grid in
%              m at each grid node of the transducer surface. Users can omit this 
%              variable, in such cases, the toolbox will automatically 
%              determine the transducer type - spherical or flat. 
%              For a spherical transducer, zGrid nodes will be set 
%              according to the sphere equation: 
%              zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xGrid.^2 - yGrid.^2) 
%              For a flat transducer, all zGrid nodes will be set to zero: 
%              zGrid = zeros(size(xGrid)).
%     'dx': x-step of the transducer Cartesian grid in m
%     'dy': y-step of the transducer Cartesian grid in m
%     'time': vector with the time samples in s at each time grid nodes 
%             for the recorded waveforms
%     'velocity': 2D or 3D matrix with the vibrational velocity 
%                          waveforms in m/s at corresponding grid nodes of
%                          the transducer
%                          if xGrid/yGrid is a 2D matrix, then 
%                          size(velocity) = [size(xGrid, 1) 
%                          size(xGrid, 2) length(time)]
%                          if xGrid/yGrid is a vector, then 
%                          size(velocity) = [length(xGrid) 
%                          length(time)]
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

%II. INPUT PARAMETERS
inputParametersFileName = '..\data_for_examples\simulation\backprojected_spherical_transient.mat'; %input with the transducer TransducerTr data structure in '.mat' format

%III. RECTANGULAR WINDOW FOR FIELD SIMULATION
xFieldBegin =  0e-3; %x, y, and z limits in m for the rectangular simulation region. Set "Begin" equal to "End" for specific coordinates to reduce the number of dimensions.
xFieldEnd   =  0e-3;
yFieldBegin = -5e-3;
yFieldEnd   =  5e-3;
zFieldBegin = 80e-3;
zFieldEnd   = 100e-3;

dxField = []; %x, y, and z grid step in m for the rectangular simulation region. Steps for singleton dimensions (where "Begin" = "End") are ignored and can be omitted.
dyField = 0.2e-3;
dzField = 0.2e-3;

%IV. TIME DOMAIN PARAMETERS
timeMax = 'auto';  % the maximum limit of the time window, which can be set either as a numeric value in s or as 'auto' to be determined automatically based on the simulation geometry.
zeroPadData = true; % perform zero-padding of the data if the edge of the input time window is less than timeMax.

%V. FREQUENCY DOMAIN PARAMETERS
powerSignificanceLevel = 0.01; %minimum significant power level of a spectral component relative to the maximum power
minFrequencySample = 'auto'; %frequency samples ranging (minFrequencySample:maxFrequencySample)*frequencyStep, where frequencyStep = 1 / "time window for TransducerTr.time"
maxFrequencySample = 'auto'; %set 'auto' to determine the values automatically based on the powerSignificanceLevel

%VI. OUTPUT VIDEO PARAMETERS
saveVideo = false; %set true to save a video with simulation results
outputVideoName = 'out_p_field.mp4'; %output filename
downsampleCoeffTime = 1; %downsampling factor for the output field related to the input 'TransducerTr.time' samples
frameRateVideo = 25; %frame rate in the output video

%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; % cascade shift for each new figure in pixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

addpath(genpath(libraryDir));

%Load input data
load(inputParametersFileName);
[TransducerTr] = check_transducer_z(TransducerTr);
[TransducerTr] = reshape_transducer_dim(TransducerTr);


%Check errors
if abs(TransducerTr.time(1))>eps
    error('The first time sample, TransducerTr.time(1), must be equal to zero!');
end
if max(diff(diff(TransducerTr.time))) > eps
    error('The time step for TransducerTr.time is nonuniform! Please check the time step!');
end

%Perform fft to extract spectral components for signals at all hologram points
expSign = 1;
timeStep = TransducerTr.time(2) - TransducerTr.time(1);
timeRinging = TransducerTr.time(end);

%Set time and frequency domain parameters
if strcmpi(timeMax,'auto')
   distMax  = sqrt((TransducerTr.xGrid(:)-xFieldEnd).^2 + (TransducerTr.yGrid(:)-yFieldEnd).^2 + (TransducerTr.zGrid(:)-zFieldEnd).^2);
   timeMax = max(distMax(:))/Medium.soundSpeed + timeRinging;
end

if zeroPadData && (exist('timeMax','var') ~= 0)
    TransducerTr.time = 0:timeStep:timeMax;
    TransducerTr.velocity = cat(3,TransducerTr.velocity, zeros([size(TransducerTr.xGrid) length(TransducerTr.time)-size(TransducerTr.velocity,3)]));
end
numberOfTimePoints = size(TransducerTr.velocity,3);
frequencyStep = (timeStep*numberOfTimePoints)^-1;

%Renormalize “MATLAB/Octave fft” into “single-sided fft”
[fArray, vSpectrum] = get_single_sided_spectrum(TransducerTr.time, TransducerTr.velocity);

%Set the automatic number of frequency samples if needed

if strcmpi(minFrequencySample,'auto')
   [minFrequencySample] = find_lowest_significant_freq(vSpectrum, powerSignificanceLevel);
end

if strcmpi(maxFrequencySample,'auto')
   maxFrequency = Medium.soundSpeed/2/max(TransducerTr.dx,TransducerTr.dy);
   maxNumFrequencySamples = round(maxFrequency/frequencyStep);
   [maxFrequencySample] = find_highest_significant_freq(frequencyStep, vSpectrum, powerSignificanceLevel, maxNumFrequencySamples, maxFrequency);
end

vSpectrum = vSpectrum(:,:,2:(maxFrequencySample+1));

if minFrequencySample > 1
    vSpectrum(:,:,1:(minFrequencySample-1)) = 0;
elseif minFrequencySample <= 0
    error('minFrequencySample must be greater than zero!');
end

%Set the radiusOfCurvature if exists
if isfield(TransducerTr, 'radiusOfCurvature')
    if ~isempty(TransducerTr.radiusOfCurvature)
        radiusOfCurvature = TransducerTr.radiusOfCurvature;
    end
end
isSphericalSource = exist('radiusOfCurvature', 'var') ~= 0;
if ~isSphericalSource
    radiusOfCurvature = [];
end

%Generate grid for simulation output
xSource = TransducerTr.xGrid;
ySource = TransducerTr.yGrid;
zSource = TransducerTr.zGrid;

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

[xField3D, yField3D, zField3D] = meshgrid(xFieldVector, yFieldVector, zFieldVector);

withinTheActiveSurface = (abs(max(TransducerTr.velocity,[],3))>eps); %surface mask of non-zero velocity
vSpectrumActiveSurface = reshape(vSpectrum, [size(vSpectrum,1)*size(vSpectrum,2), size(vSpectrum,3)]);
vSpectrumActiveSurface = vSpectrumActiveSurface(withinTheActiveSurface,:);

SourceParameters = [];
SourceParameters.xGrid = xSource(withinTheActiveSurface);
SourceParameters.yGrid = ySource(withinTheActiveSurface);
SourceParameters.zGrid = zSource(withinTheActiveSurface);
SourceParameters.dx = TransducerTr.dx;
SourceParameters.dy = TransducerTr.dy;
SourceParameters.input = vSpectrumActiveSurface;

FieldParameters = [];
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;

isTransient = true; %set to true for the transient regime, false for the single-frequency regime

if isSphericalSource
    regime = 3; % 3 Forward-projection: V on a sphere --> P at an arbitrary set of points
    %Rayleigh simulator function with the complex pressure amplitude output
    [ pSpectrumField ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

else

    regime = 4; % 4 Forward-projection: V on a plane --> P at an arbitrary set of points
    %Rayleigh simulator function with the complex pressure amplitude output
    [ pSpectrumField ] = rayleigh_simulator(expSign, frequencyStep, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

end

pSpectrumField = reshape(pSpectrumField, [length(yFieldVector) length(xFieldVector) length(zFieldVector) maxFrequencySample]);

%Renormalize single-sided fft into MATLAB/Octave fft
pSpectrumFieldMatlabCalibration = pSpectrumField*numberOfTimePoints;
clear pSpectrumField;

if mod(numberOfTimePoints,2) == 0
    pSpectrumFieldMatlabCalibration(:,:,:,2:(end-1)) = pSpectrumFieldMatlabCalibration(:,:,:,2:(end-1))/2;
else
    pSpectrumFieldMatlabCalibration(:,:,:,2:end) = pSpectrumFieldMatlabCalibration(:,:,:,2:end)/2;
end

%Zero-pad the obtained spectrum to maintain the same time step of the output as in input time grid 'TransducerTr.time'
sizeOutput = size(pSpectrumFieldMatlabCalibration);
lastDimOutput = length(sizeOutput);
pSpectrumFieldTotal = cat(lastDimOutput, zeros(sizeOutput(1:end-1)),...
    pSpectrumFieldMatlabCalibration,...
    zeros([sizeOutput(1:end-1) (numberOfTimePoints-2*maxFrequencySample-1)]),...
    flip(conj(pSpectrumFieldMatlabCalibration),lastDimOutput));
clear pSpectrumFieldMatlabCalibration;

%Back to time domain
pFieldFrames3D = real(ifft(pSpectrumFieldTotal,[],lastDimOutput));
clear pSpectrumFieldTotal;

distMin  = min(sqrt((TransducerTr.xGrid(:)-xFieldBegin).^2 + (TransducerTr.yGrid(:)-yFieldBegin).^2 + (TransducerTr.zGrid(:)-zFieldBegin).^2));
distMax  = max(sqrt((TransducerTr.xGrid(:)-xFieldEnd).^2 + (TransducerTr.yGrid(:)-yFieldEnd).^2 + (TransducerTr.zGrid(:)-zFieldEnd).^2));

iBegin = round((distMin/Medium.soundSpeed - timeRinging)/timeStep);
iEnd =   round((distMax/Medium.soundSpeed + timeRinging)/timeStep);

if iBegin < 1
  iBegin = 1;
end

if iEnd > length(TransducerTr.time)
    iEnd = length(TransducerTr.time);
end

%GUI part to show results in a slice-by-slice format and save them
if abs(ndims(squeeze(pFieldFrames3D)) - 3) > eps
    warning('The output field is not 2D, so it can''t be displayed in a frame-by-frame format. Users can access the data using the workspace variables xField3D, yField3D, zField3D, and pFieldFrames3D and display them as desired.');
else
    [outX, outY, pField2D, outputAxes, constAxis, constAxisValue] =  gui_fp_transient_frame_by_frame(xField3D,yField3D, zField3D, pFieldFrames3D, timeStep, iBegin, iEnd, true, radiusOfCurvature);


    if saveVideo
        figure;
        v = VideoWriter(outputVideoName, 'MPEG-4');
        set(v,'FrameRate',frameRateVideo);
        set(v,'Quality',100);
        open(v);
        pLimits = [min(real(pFieldFrames3D(:))) max(pFieldFrames3D(:))];

        for i = iBegin:downsampleCoeffTime:iEnd

            pFrameCurrent = squeeze(pField2D(:,:,i));

            imagesc(outX([1 end])*1e3, outY([1 end])*1e3, pFrameCurrent);
            set(gca,'YDir','normal');

            axis equal

            caxis(pLimits);

            colormap jet
            h = colorbar;
            axis tight;

            title(['Pressure ' outputAxes '-distribution in Pa (' constAxis ' = ' num2str(constAxisValue*1e3) ' mm) at {\itt} =' num2str(TransducerTr.time(i)*1e6, '%1.4f') ' us' ]);

            xlabel('z, mm');
            ylabel('x, mm');
            frame_out = getframe(gcf);
            writeVideo(v,frame_out);

        end
        close(v);
    end
end

disp('xDDx simulation completed!');