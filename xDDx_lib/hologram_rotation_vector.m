% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xMax, yMax, zMax, directionVector, RotationLine] = hologram_rotation_vector(xField3D, yField3D, zField3D, fieldOut3D, izMaxArray)

fieldOutReshaped = reshape(fieldOut3D(:,:,izMaxArray),[size(fieldOut3D,1)*size(fieldOut3D,2),length(izMaxArray)]);
[~,iMaxLayers] = max(abs(fieldOutReshaped),[],1);
x = reshape(xField3D(:,:,izMaxArray),size(fieldOutReshaped));
y = reshape(yField3D(:,:,izMaxArray),size(fieldOutReshaped));
z = reshape(zField3D(:,:,izMaxArray),size(fieldOutReshaped));

[~,iMax] = max(abs(fieldOut3D(:)));
xMax = xField3D(iMax);
yMax = yField3D(iMax);
zMax = zField3D(iMax);

xCenters = zeros(1,length(izMaxArray));
yCenters = zeros(1,length(izMaxArray));
zCenters = zeros(1,length(izMaxArray));
for j = 1:size(fieldOutReshaped,2)
    xCenters(j) = x(iMaxLayers(j), j);
    yCenters(j) = y(iMaxLayers(j), j);
    zCenters(j) = z(iMaxLayers(j), j);
end
xyzCenters = [xCenters; yCenters; zCenters];

xyzMean = mean(xyzCenters,2);
A = (xyzCenters-xyzMean);

[V,D] = eig((A').'*(A'),'vector');
[~,iMaxLayers] = max(D);
directionVector = V(:,iMaxLayers);

if directionVector(3) < 0
 directionVector = -directionVector;
end

RotationLine = [];
RotationLine.angleZ = acosd(dot(directionVector,[0; 0; 1]));

t = directionVector'*A;
tBegin = min(t);
tEnd = max(t);
RotationLine.xyzLine = xyzMean + [tBegin, tEnd].*directionVector; % size 3x2

RotationLine.xyzInputPoints = xyzCenters;

end