function MaterialMatrixSimulation = generate_simulation_model(radialReserveX, radialReserveY, aperture, radiusOfCurvature, ...
    zGridMediumVec, izBoundaryCondition, dx, dy, dz, izTarget, zFieldEnd, ...
    xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, ~, ...
    pml_x_size, pml_y_size, pml_z_size, ixTarget, iyTarget, ...
    nxMedium, nyMedium, nzMedium, MaterialMatrix, ...
    soundSpeedContactMedium, densityContactMedium, alphaContactMedium, strongMemorySavingMode)
%GENERATE_SIMULATION_MODEL Generate simulation grid and material matrix
%   MaterialMatrixSimulation = generate_simulation_model(...)
%   Calculates the simulation grid dimensions, padding parameters, and 
%   resized material matrix for kWave simulation.
%
%   Inputs:
%       radialReserveX, radialReserveY - Reserve factors for the boundary
%                  condition size in the x and y directions
%       aperture - Nominal transducer aperture in m. Pass a scalar for the
%                  same x/y aperture, or [apertureX apertureY].
%       radiusOfCurvature - Transducer radius of curvature (m)
%       zGridMediumVec - Z grid vector for medium (m)
%       izBoundaryCondition - Index of boundary condition in z dimension
%       dx, dy, dz - Grid spacing in each dimension (m)
%       izTarget - Index of target in z dimension
%       zFieldEnd - End position of field in z dimension (m)
%       xFieldBegin, xFieldEnd - x-limits of requested output window (m)
%       yFieldBegin, yFieldEnd - y-limits of requested output window (m)
%       zFieldBegin - z-begin of requested output window (m)
%       pml_x_size, pml_y_size, pml_z_size - PML sizes in each dimension
%       ixTarget, iyTarget - Target indices in x and y dimensions
%       nxMedium, nyMedium, nzMedium - Medium grid dimensions
%       MaterialMatrix - Material matrix structure
%       soundSpeedContactMedium - Sound speed of contact medium (m/s)
%       densityContactMedium - Density of contact medium (kg/m^3)
%       alphaContactMedium - Attenuation of contact medium (dB/cm/MHz)
%
%   Output:
%       MaterialMatrixSimulation - Material matrix structure defined for the simulation grid

    if nargin < 29
        strongMemorySavingMode = false;
    end

    isSphericalSource = ~isempty(radiusOfCurvature);
    [apertureX, apertureY, apertureScalar] = parse_aperture(aperture);

    if isSphericalSource
        transverseSizeBoundaryConditionX = (1 + radialReserveX) * apertureX * ...
            (radiusOfCurvature - zGridMediumVec(izBoundaryCondition)) / ...
            sqrt(radiusOfCurvature^2 - (apertureScalar/2)^2);
        transverseSizeBoundaryConditionY = (1 + radialReserveY) * apertureY * ...
            (radiusOfCurvature - zGridMediumVec(izBoundaryCondition)) / ...
            sqrt(radiusOfCurvature^2 - (apertureScalar/2)^2);

        nxSimulationGridTemp = round(transverseSizeBoundaryConditionX / dx);
        nySimulationGridTemp = round(transverseSizeBoundaryConditionY / dy);
        nzSimulationGridTemp = izTarget + round((zFieldEnd - radiusOfCurvature) / dz) - ...
            izBoundaryCondition + 1;
    else
        % Flat transducer: source is at izBoundaryCondition and centered at target.
        transverseSizeBoundaryConditionX = (1 + radialReserveX) * apertureX;
        transverseSizeBoundaryConditionY = (1 + radialReserveY) * apertureY;

        % Expand by requested 3D output window when needed.
        windowSizeX = abs(xFieldEnd - xFieldBegin);
        windowSizeY = abs(yFieldEnd - yFieldBegin);
        transverseSizeBoundaryConditionX = max(transverseSizeBoundaryConditionX, windowSizeX);
        transverseSizeBoundaryConditionY = max(transverseSizeBoundaryConditionY, windowSizeY);

        nxSimulationGridTemp = round(transverseSizeBoundaryConditionX / dx);
        nySimulationGridTemp = round(transverseSizeBoundaryConditionY / dy);

        zBoundary = zGridMediumVec(izBoundaryCondition);
        nzSimulationGridTemp = round((zFieldEnd - zBoundary) / dz) + 1;
    end
    
    % Adjust grid sizes for optimal FFT performance
    [nxSimulationGrid, nySimulationGrid, nzSimulationGrid] = adjust_fft_size(...
        nxSimulationGridTemp, nySimulationGridTemp, nzSimulationGridTemp,...
        pml_x_size, pml_y_size, pml_z_size);
    
    % Calculate simulation grid indices
    ixBeginSimulationGrid = ixTarget - (nxSimulationGrid / 2);
    ixEndSimulationGrid = ixTarget + (nxSimulationGrid / 2 - 1);
    iyBeginSimulationGrid = iyTarget - (nySimulationGrid / 2);
    iyEndSimulationGrid = iyTarget + (nySimulationGrid / 2 - 1);
    izBeginSimulationGrid = izBoundaryCondition;
    izEndSimulationGrid = izBoundaryCondition + (nzSimulationGrid - 1);
    
    % Calculate padding requirements
    nxPadBegin = (1 - ixBeginSimulationGrid);
    nyPadBegin = (1 - iyBeginSimulationGrid);
    nzPadBegin = (1 - izBeginSimulationGrid);
    
    nxPadEnd = (ixEndSimulationGrid - nxMedium);
    nyPadEnd = (iyEndSimulationGrid - nyMedium);
    nzPadEnd = (izEndSimulationGrid - nzMedium);
    
    % Resize material matrix
    if strongMemorySavingMode
        MaterialMatrixSimulation = resize_material_matrix_strong(MaterialMatrix, ...
            nxPadBegin, nyPadBegin, nzPadBegin, ...
            nxPadEnd, nyPadEnd, nzPadEnd, ...
            soundSpeedContactMedium, densityContactMedium, alphaContactMedium);
    else
        MaterialMatrixSimulation = resize_material_matrix(MaterialMatrix, ...
            nxPadBegin, nyPadBegin, nzPadBegin, ...
            nxPadEnd, nyPadEnd, nzPadEnd, ...
            soundSpeedContactMedium, densityContactMedium, alphaContactMedium);
    end
    
end

function [apertureX, apertureY, apertureScalar] = parse_aperture(aperture)
    if ~isnumeric(aperture) || isempty(aperture) || any(aperture(:) <= 0)
        error('aperture must be a positive scalar or a two-element positive numeric vector.');
    end
    if isscalar(aperture)
        apertureX = aperture;
        apertureY = aperture;
    elseif numel(aperture) == 2
        apertureX = aperture(1);
        apertureY = aperture(2);
    else
        error('aperture must be a positive scalar or a two-element positive numeric vector.');
    end
    apertureScalar = max(apertureX, apertureY);
end

