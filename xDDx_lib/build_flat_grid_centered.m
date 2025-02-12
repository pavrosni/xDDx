% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
% 
function [xFlatBoundary, yFlatBoundary] = build_flat_grid_centered(nxFlatBoundary,dxFlatBoundary,nyFlatBoundary,dyFlatBoundary)

if mod(nxFlatBoundary,2) ~= 0
    xFlatBoundaryVector = ((-(nxFlatBoundary-1)/2):((nxFlatBoundary-1)/2)) * dxFlatBoundary;
else
    xFlatBoundaryVector = ((-nxFlatBoundary/2):(nxFlatBoundary/2-1)) * dxFlatBoundary;
end

if mod(nyFlatBoundary,2) ~= 0
    yFlatBoundaryVector = ((-(nyFlatBoundary-1)/2):((nyFlatBoundary-1)/2)) * dyFlatBoundary;
else
    yFlatBoundaryVector = ((-nyFlatBoundary/2):(nyFlatBoundary/2-1)) * dyFlatBoundary;
end

[xFlatBoundary, yFlatBoundary] = meshgrid(xFlatBoundaryVector, yFlatBoundaryVector);

end

