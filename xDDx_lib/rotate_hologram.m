% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xHoloRotated, yHoloRotated, zHoloRotated, holoRotated] = rotate_hologram(xMax, yMax, zMax, directionVector, expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, MediumParameters, ServiceParameters, zMaxAcoustCoord, shiftHoloStep)

xHolo = SourceParameters.xGrid;
yHolo = SourceParameters.yGrid;
[xHoloGlobal, yHoloGlobal, zHoloGlobal, zPositionHoloRotated] = calculate_rotated_coordinates(xMax, yMax, zMax, directionVector, xHolo, yHolo, zMaxAcoustCoord, shiftHoloStep);

FieldParameters.xGrid = xHoloGlobal;
FieldParameters.yGrid = yHoloGlobal;
FieldParameters.zGrid = zHoloGlobal;
[ holoRotated ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, MediumParameters, ServiceParameters);

xHoloRotated = xHolo;
yHoloRotated = yHolo;
zHoloRotated = zPositionHoloRotated*ones(size(xHoloRotated));

end

