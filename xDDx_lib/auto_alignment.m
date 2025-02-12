% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, RotationLine, fvFocalLobeCoarse, fvFocalLobeFine] = auto_alignment(simulationDevice, Geometry, HologramSf, Medium, ServiceParameters, focalLobeLevel, ppwRadialCoarse3D, ppwAxialCoarse3D, ppwRadialFine3D, ppwAxialFine3D, shiftHoloStep)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
expSign = HologramSf.expSign;
frequency = HologramSf.frequency;

%Approximating the size of the focal lobe given the geometry of the transducer
%For equations see Rosnitskiy et al., IEEE UFFC, 64(2), PP. 374 - 390, 2018
%DOI: 10.1109/TUFFC.2016.2619913
radiusOfCurvature = Geometry.radiusOfCurvature;
aperture = Geometry.apertureMin;

xHolo = HologramSf.xGrid;
yHolo = HologramSf.yGrid;
zHolo = zeros(size(xHolo));

G = (2*pi*frequency/Medium.soundSpeed)*(aperture/2)^2/2/radiusOfCurvature;
rBesselRoot = fzero(@(rBessel) besselj(1, rBessel), 3);
rRoot1 = -2*rBesselRoot*aperture/4/G;
rRoot2 =  2*rBesselRoot*aperture/4/G;
funzRoot1 = @(fNumber,k,F) 4*pi*fNumber/k * (k*F*(-k*F-pi)*sqrt(4*fNumber^2-1)-2*fNumber*((k*F)^2 + 3*pi*k*F + 2*pi^2)) / ((k*F)^2 + 16*pi*(pi + k*F)*fNumber^2) + F;
funzRoot2 = @(fNumber,k,F) 4*pi*fNumber/k * (k*F*(+k*F-pi)*sqrt(4*fNumber^2-1)+2*fNumber*((k*F)^2 - 3*pi*k*F + 2*pi^2)) / ((k*F)^2 + 16*pi*(pi - k*F)*fNumber^2) + F;

zRoot1 = funzRoot1(radiusOfCurvature/aperture, 2*pi*frequency/Medium.soundSpeed, radiusOfCurvature);
zRoot2 = funzRoot2(radiusOfCurvature/aperture, 2*pi*frequency/Medium.soundSpeed, radiusOfCurvature);

zPrefocalShiftApprox = radiusOfCurvature - HologramSf.zPosition;
xFieldMin =  rRoot1-HologramSf.errorPositioning;
xFieldMax =  rRoot2+HologramSf.errorPositioning;
yFieldMin =  rRoot1-HologramSf.errorPositioning;
yFieldMax =  rRoot2+HologramSf.errorPositioning;
zFieldMin = zPrefocalShiftApprox - (radiusOfCurvature - zRoot1) - HologramSf.errorPositioning;
zFieldMax = zPrefocalShiftApprox + (zRoot2 - radiusOfCurvature) + HologramSf.errorPositioning;

if zFieldMin < 0
    zFieldMin = shiftHoloStep;
end

% Creating a coarse 3D grid for the "First Iteration" focal lobe simulation
if ~isOctave
    dxFieldCoarse = round(Medium.soundSpeed/frequency/ppwRadialCoarse3D,1,'significant');
    dyFieldCoarse = round(Medium.soundSpeed/frequency/ppwRadialCoarse3D,1,'significant');
    dzFieldCoarse = round(Medium.soundSpeed/frequency/ppwAxialCoarse3D,1,'significant');
else
    dxFieldCoarse = Medium.soundSpeed/frequency/ppwRadialCoarse3D;
    dyFieldCoarse = Medium.soundSpeed/frequency/ppwRadialCoarse3D;
    dzFieldCoarse = Medium.soundSpeed/frequency/ppwAxialCoarse3D;
end

[xField3D, yField3D, zField3D] = meshgrid(xFieldMin:dxFieldCoarse:xFieldMax,...
    yFieldMin:dyFieldCoarse:yFieldMax,...
    zFieldMin:dzFieldCoarse:zFieldMax);

%"First Iteration" coarse-grid focal lobe simulation
SourceParameters.xGrid = xHolo;
SourceParameters.yGrid = yHolo;
SourceParameters.zGrid = zHolo;
SourceParameters.dx = HologramSf.dx;
SourceParameters.dy = HologramSf.dy;
SourceParameters.input = HologramSf.complexPressureAmplitude;
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;
regime = 5;

isTransient = false; %set to true for the transient regime, false for the single-frequency regime

[ pOutCoarse3D ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);


%Creating an isosurface for the coarse-grid focal lobe
focalLobeLevelCoarse = 0.9*focalLobeLevel; % 10 percent of isovalue is reserved for the fine-grid simulation
fvFocalLobeCoarse = isosurface(xField3D,yField3D,zField3D,abs(pOutCoarse3D)/max(abs(pOutCoarse3D(:))), focalLobeLevelCoarse);

%Plot the coarse-grid focal lobe isosurface
positionCurrent = get(groot,'DefaultFigurePosition');

%Creating a fine 3D grid for the "Second Iteration" focal lobe simulation
xFieldMin = min(fvFocalLobeCoarse.vertices(:,1));
xFieldMax = max(fvFocalLobeCoarse.vertices(:,1));
yFieldMin = min(fvFocalLobeCoarse.vertices(:,2));
yFieldMax = max(fvFocalLobeCoarse.vertices(:,2));
zFieldMin = min(fvFocalLobeCoarse.vertices(:,3));
zFieldMax = max(fvFocalLobeCoarse.vertices(:,3));

if ~isOctave
    dxField = round(Medium.soundSpeed/frequency/ppwRadialFine3D,1,'significant');
    dyField = round(Medium.soundSpeed/frequency/ppwRadialFine3D,1,'significant');
    dzField = round(Medium.soundSpeed/frequency/ppwAxialFine3D,1,'significant');
else
    dxField = Medium.soundSpeed/frequency/ppwRadialFine3D;
    dyField = Medium.soundSpeed/frequency/ppwRadialFine3D;
    dzField = Medium.soundSpeed/frequency/ppwAxialFine3D;
end

[xField3D, yField3D, zField3D] = meshgrid(xFieldMin:dxField:xFieldMax,...
    yFieldMin:dyField:yFieldMax,...
    zFieldMin:dzField:zFieldMax);

%"Second Iteration" fine-grid focal lobe simulation
SourceParameters.xGrid = xHolo;
SourceParameters.yGrid = yHolo;
SourceParameters.zGrid = zHolo;
SourceParameters.dx = HologramSf.dx;
SourceParameters.dy = HologramSf.dy;
SourceParameters.input = HologramSf.complexPressureAmplitude;
FieldParameters.xGrid = xField3D;
FieldParameters.yGrid = yField3D;
FieldParameters.zGrid = zField3D;
regime = 5;
[ pOut3D ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

%Get the focal lobe isosurface
fvFocalLobeFine = isosurface(xField3D,yField3D,zField3D,abs(pOut3D)/max(abs(pOut3D(:))), focalLobeLevel);

%Using 'fvFocalLobeFine' to find rotation angle of the hologram relative to the transducer's axis of symmetry
zForRotation1 = min(fvFocalLobeFine.vertices(:,3));
zForRotation2 = max(fvFocalLobeFine.vertices(:,3));
iz1 = find(zField3D(1,1,:) >= zForRotation1, 1, 'first');
iz2 = find(zField3D(1,1,:) <= zForRotation2, 1, 'last');
pOut3D = pOut3D(:,:,iz1:iz2);
xField3D = xField3D(:,:,iz1:iz2);
yField3D = yField3D(:,:,iz1:iz2);
zField3D = zField3D(:,:,iz1:iz2);
layersToProcess3D = 1:size(pOut3D,3);
[xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, RotationLine] = hologram_rotation_vector(xField3D, yField3D, zField3D, pOut3D, layersToProcess3D);

end