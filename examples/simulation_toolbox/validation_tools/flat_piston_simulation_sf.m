% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% 
% This script demonstrates simulation precision using an example of an 
% idealized piston flat source (uniform distribution of the vibrational 
% velocity over the surface). The 3D field of the source is simulated 
% numerically using the toolbox, and the numerical result along the 
% transducer axis is compared to the analytical solution 
% [O'Neil, DOI: 10.1121/1.1906542].

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


%III. TRANSDUCER PARAMETERS
aperture = 20e-3; %in m
frequency = 1e6; %in Hz
sourceStepX = 0.7e-3;  %x grid step for source representation in m
sourceStepY = 0.7e-3;  %y grid step for source representation in m
initialVelocity = 1/Medium.soundSpeed/Medium.density; %vibrational velocity at the surface of the source in m/s


%IV. 3D WINDOW FOR FIELD SIMULATION
xFieldBegin = -7e-3; %x, y, and z limits in m for the rectangular 3D simulation region
xFieldEnd   =  7e-3;
yFieldBegin = -7e-3;
yFieldEnd   =  7e-3;
zFieldBegin = 10e-3;
zFieldEnd   = 120e-3;

dxField = 0.25e-3; %x, y, and z grid step in m for the rectangular 3D simulation region
dyField = 0.25e-3;
dzField = 0.5e-3;


%V. ISOLEVEL PARAMETERS
levelArray = [0.8 0.5 0.3 0.1]; %isolevels related to the pressure maximum for extracting the isosurface of the simulated 3D field
transparencyArray = [0.6 0.5 0.2 0.05]; %transparency of the isosurfaces for the levels from levelArray


%SERVICE PARAMETERS
ServiceParameters.threadsPerBlockGPU = 128; %number of threads per block for GPU (if applicable)
figuresCascadeShift = 50; %cascade shift for each new figure in pixels
colorActiveSurface = [1 0 0]; %RGB color of the active surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

addpath(genpath(libraryDir));


%Generate the Transducer grid that includes the center of symmetry of the transducer

%Reserve one grid points from each boundary of the source grid
nxSource = round(aperture/sourceStepX)+2;
nySource = round(aperture/sourceStepY)+2;

%Build the boundary condition grid that includes the center of symmetry of the transducer
[xSource, ySource] = build_flat_grid_centered(nxSource,sourceStepX,nySource,sourceStepY);
zSource = zeros(size(xSource));

complexVelocityAmplitude = initialVelocity*ones(size(xSource));
complexVelocityAmplitude(xSource.^2 + ySource.^2 > (aperture/2)^2) = 0;

Transducer = [];
Transducer.expSign = 1;
Transducer.frequency = frequency;
Transducer.xGrid = xSource;
Transducer.yGrid = ySource;
Transducer.zGrid = zSource;
Transducer.dx = sourceStepX;
Transducer.dy = sourceStepY;
Transducer.complexVelocityAmplitude = complexVelocityAmplitude;

%Generate 3D grid for simulation output
xFieldVector = xFieldBegin : dxField : xFieldEnd;
yFieldVector = yFieldBegin : dyField : yFieldEnd;
zFieldVector = zFieldBegin : dzField : zFieldEnd;

[xField3D, yField3D, zField3D] = meshgrid(xFieldVector, yFieldVector, zFieldVector);

%Forwardpropagate the hologram to calculate the complex pressure amplitude at the nodes of the grid
expSign = Transducer.expSign;

withinTheActiveSurface = (abs(Transducer.complexVelocityAmplitude)>eps); %surface mask of non-zero velocity 

SourceParameters = [];
SourceParameters.xGrid = Transducer.xGrid(withinTheActiveSurface);
SourceParameters.yGrid = Transducer.yGrid(withinTheActiveSurface);
SourceParameters.zGrid = Transducer.zGrid(withinTheActiveSurface);
SourceParameters.dx = Transducer.dx;
SourceParameters.dy = Transducer.dy;
SourceParameters.input = Transducer.complexVelocityAmplitude(withinTheActiveSurface);

FieldParameters = [];
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;

regime = 4; % 4 Forwardpropagation: V on a plane --> P at an arbitrary set of points

isTransient = false;

%Rayleigh simulator function with the complex pressure amplitude output
[ pField3D ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

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
imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,abs(complexVelocityAmplitude)*1e3);
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
imagesc(xSource([1 end])*1e3,ySource([1 end])*1e3,angle(complexVelocityAmplitude));
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


%Plot 3D field and source surface
figure;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
positionCurrent = positionCurrent + [figuresCascadeShift -figuresCascadeShift 0 0];
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
titleFor3D = ['Pressure amplitude max is ' num2str(pMax3D) ' Pa at ' '(' num2str(xMax3D*1e3) ', ' num2str(yMax3D*1e3) ', ' num2str(zMax3D*1e3) ') mm' ];
if ~isOctave
    cbh = colorbar;
    caxis([0 1]);
    set(cbh,'XTick',flip(levelArray));
else
     titleFor3D = {titleFor3D, ['pressure amplitude iso-levels: ' num2str(levelArray)]};
end
title(titleFor3D);


notActiveSurface = abs(withinTheActiveSurface-1);
activeSurfaceImage = zeros(size(withinTheActiveSurface,1),size(withinTheActiveSurface,2),3);
activeSurfaceImage(:,:,1) = colorActiveSurface(1)*withinTheActiveSurface;
activeSurfaceImage(:,:,2) = colorActiveSurface(2)*withinTheActiveSurface;
activeSurfaceImage(:,:,3) = colorActiveSurface(3)*withinTheActiveSurface;
activeSurfaceImage = activeSurfaceImage + notActiveSurface;
zSourceToPlot = zSource;
zSourceToPlot(~withinTheActiveSurface) = NaN;
surface(xSource*1e3, ySource*1e3, zSourceToPlot*1e3, 'FaceColor', 'r','EdgeColor','none');
axis equal;
if ~isOctave
    camlight;
else
    zlim([0 max(zField3D(:))*1e3]);
end
view(3);

%GUI part to show 3D results in a slice-by-slice format and save them
gui_plot_3d(xField3D,yField3D, zField3D, pField3D, levelArray, transparencyArray, false);


%Calculate the field at the axis of the source using the analytical O'Neil's solution
z = zFieldVector;
p0 = Medium.density*Medium.soundSpeed*initialVelocity;
a = aperture/2;
k = 2*pi*frequency/Medium.soundSpeed;
pFieldAnalytical = 2*1i*p0*exp(-1i*k/2*(sqrt(a^2+z.^2) + z)).*sin(k/2*(sqrt(a^2+z.^2) - z));

%Extract the field at the axis of the source from the numerical simulation results
ixZero = find(abs(xFieldVector) < eps,1,'first');
iyZero = find(abs(yFieldVector) < eps,1,'first');
pFieldNumerical = squeeze(pField3D(iyZero,ixZero,:));

%Compare analytical and numerical solutions
figure;
hold on;
set(gcf, 'Units', 'pixels', 'Position', positionCurrent);
plot(z*1e3, abs(pFieldNumerical),'LineWidth',2);
plot(z*1e3, abs(pFieldAnalytical),'--','LineWidth',2);
xlim([z(1) z(end)]*1e3);
xlabel('z, mm');
ylabel('Pressure Amplitude, Pa');
title(['Analytical Solution vs. Numerical with spatial steps of ' num2str(sourceStepX*1e3) ' by ' num2str(sourceStepY*1e3) ' mm']);
legend('Numerical', 'Analytical');
grid on;
grid minor;

disp('xDDx simulation completed!');