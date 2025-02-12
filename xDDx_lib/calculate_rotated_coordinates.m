% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xHoloRotatedMech, yHoloRotatedMech, zHoloRotatedMech, zPositionHoloRotatedAcoust] = calculate_rotated_coordinates(xMax, yMax, zMax, directionVector, xHolo, yHolo, zMaxAcoust, shiftHoloStep)
ort1Mech = [1; 0; 0];
ort2Mech = [0; 1; 0];
ort3Mech = [0; 0; 1];
ort1Acoust = [(1+directionVector(1)^2/directionVector(3)^2)^-0.5; 0; -directionVector(1)/directionVector(3)*(1+directionVector(1)^2/directionVector(3)^2)^-0.5];
if ort1Acoust(1) < 0
    ort1Acoust = -ort1Acoust;
end
ort1Acoust = ort1Acoust/norm(ort1Acoust);
ort3Acoust = directionVector;
ort2Acoust = cross(ort3Acoust, ort1Acoust);
ort2Acoust = ort2Acoust/norm(ort2Acoust);
needToShiftHolo = true;
shiftRotated = 0;

while needToShiftHolo
    shiftRotated = shiftRotated + shiftHoloStep;
    shiftFromFocusForRotated = zMaxAcoust-shiftRotated;
    
    xCenterAtNewHoloMech = xMax - shiftFromFocusForRotated*directionVector(1);
    yCenterAtNewHoloMech = yMax - shiftFromFocusForRotated*directionVector(2);
    zCenterAtNewHoloMech = zMax - shiftFromFocusForRotated*directionVector(3);
    
    zHolo = zeros(size(xHolo));
    
    xHoloRotatedMech = xCenterAtNewHoloMech*ort1Mech(1) + yCenterAtNewHoloMech*ort2Mech(1) + zCenterAtNewHoloMech*ort3Mech(1) + xHolo*ort1Acoust(1) + yHolo*ort2Acoust(1) + zHolo*ort3Acoust(1);
    yHoloRotatedMech = xCenterAtNewHoloMech*ort1Mech(2) + yCenterAtNewHoloMech*ort2Mech(2) + zCenterAtNewHoloMech*ort3Mech(2) + xHolo*ort1Acoust(2) + yHolo*ort2Acoust(2) + zHolo*ort3Acoust(2);
    zHoloRotatedMech = xCenterAtNewHoloMech*ort1Mech(3) + yCenterAtNewHoloMech*ort2Mech(3) + zCenterAtNewHoloMech*ort3Mech(3) + xHolo*ort1Acoust(3) + yHolo*ort2Acoust(3) + zHolo*ort3Acoust(3);
        
    needToShiftHolo = (min(zHoloRotatedMech(:)) <= shiftHoloStep);
end

zPositionHoloRotatedAcoust = zMaxAcoust - shiftFromFocusForRotated;

end

