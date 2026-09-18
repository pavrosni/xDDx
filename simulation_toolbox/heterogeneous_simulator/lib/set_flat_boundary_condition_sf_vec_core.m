function [vFlatBoundary] = ...
    set_flat_boundary_condition_sf_vec_core(TransducerSf, soundSpeed, density, simulationDevice, ...
                                        xFlatBoundary, yFlatBoundary, zFlatBoundary, dxFlatBoundary, dyFlatBoundary, ...
                                        shiftHoloDistInWl, ServiceParameters)
%SET_FLAT_BOUNDARY_CONDITION_SF_CORE  Build flat boundary condition for a spherical transducer.
%
%   This function back-projects the vibrational velocity defined on a
%   spherical transducer surface (given by TransducerSf) to a flat plane
%   located at the apex of the transducer (z = 0). The result is the
%   complex vibrational velocity vFlatBoundary on a Cartesian grid
%   (xFlatBoundary, yFlatBoundary, zFlatBoundary).
%
%   INPUT:
%       TransducerSf      - struct with spherical transducer parameters:
%                           .expSign
%                           .frequency
%                           .radiusOfCurvature
%                           .xGrid, .yGrid, .zGrid
%                           .dx, .dy
%                           .complexVelocityAmplitude
%       soundSpeed        - sound speed (m/s)
%       density           - density (kg/m^3)
%                           .soundSpeed, .density
%       simulationDevice  - 'cuda' or 'cpu'
%       xFlatBoundary     - x-coordinates (m) of the flat BC grid (matrix)
%       yFlatBoundary     - y-coordinates (m) of the flat BC grid (matrix)
%       shiftHoloDistInWl - (optional) minimum back‑projection distance in wavelengths
%                           (default: 10)
%       ServiceParameters - (optional) struct with technical parameters for the
%                           Rayleigh simulator (e.g., threadsPerBlockGPU) (default: ServiceParameters.threadsPerBlockGPU = 128)
%
%   OUTPUT:
%       xFlatBoundary, yFlatBoundary, zFlatBoundary - coordinates (m) of
%           the flat boundary condition grid (Cartesian)
%       vFlatBoundary - complex vibrational velocity (m/s) at the flat
%           boundary nodes
%
%   The implementation is based on the original reference script
%   reference_scripts/set_flat_boundary_condition_sf.m, but without any
%   plotting or 3D field simulation.


    if nargin < 11
      ServiceParameters = [];
      ServiceParameters.threadsPerBlockGPU = 128;
    end

    if nargin < 10
        shiftHoloDistInWl = 10;
    end

    Medium = [];
    Medium.soundSpeed = soundSpeed;
    Medium.density = density;

    % Ensure grids are 2D if a vector representation is provided
    if isvector(TransducerSf.xGrid)
        TransducerSf = reshape_transducer_dim(TransducerSf);
    end

    activeTransducer = abs(TransducerSf.complexVelocityAmplitude) > eps('single');
    if ~any(activeTransducer(:))
        error(['All transducer grid points have negligible complexVelocityAmplitude ', ...
            '(|v| <= eps(''single'')).']);
    end

    % Geometric and frequency parameters
    expSign   = TransducerSf.expSign;
    frequency = TransducerSf.frequency;

    if isfield(TransducerSf, 'radiusOfCurvature') && ~isempty(TransducerSf.radiusOfCurvature)
        radiusOfCurvature = TransducerSf.radiusOfCurvature;
    else
        error('TransducerSf.radiusOfCurvature must be provided for a spherical source.');
    end

    isSphericalSource = ~isempty(radiusOfCurvature);
    if ~isSphericalSource
        error('Your transducer is already flat. No need to apply the boundary condition transfer.');
    end

    % Derive grid sizes and steps from the provided boundary grids
    [nyFlatBoundary, nxFlatBoundary] = size(xFlatBoundary);


    % Minimum back‑projection distance in meters
    shiftHoloDist = shiftHoloDistInWl * Medium.soundSpeed / frequency;

    % Build the boundary condition grid that includes the center of symmetry
    %zFlatBoundary = zeros(size(xFlatBoundary));

    % Back‑project the source vibrational velocity to calculate the complex
    % pressure amplitude at the nodes of the flat boundary condition at z = 0.
    % If the source is too close to z = 0, additional back/forward projection
    % is performed to satisfy shiftHoloDistInWl.

    if (min(TransducerSf.zGrid(:))-zFlatBoundary(1) < shiftHoloDist)

        % Save parameters of the desired boundary condition at z = 0
        xFlatBoundaryDesired = xFlatBoundary;
        yFlatBoundaryDesired = yFlatBoundary;
        zFlatBoundaryDesired = zFlatBoundary;



        % Increase the size of the shifted flat boundary condition as it is
        % more distant from the focus than the desired one

        % Detect whether input grids are ndgrid (x along rows) or meshgrid (x along columns)
        isNdgrid = any(abs(diff(xFlatBoundary(:, 1))) > eps('single'));
        if isNdgrid
            xFlatBoundaryVector = squeeze(xFlatBoundary(:, 1));
            yFlatBoundaryVector = squeeze(yFlatBoundary(1, :));
        else
            xFlatBoundaryVector = squeeze(xFlatBoundary(1, :));
            yFlatBoundaryVector = squeeze(yFlatBoundary(:, 1));
        end

        maxTransverseSize = max(abs(xFlatBoundaryVector(end)-xFlatBoundaryVector(1)), ...
                                abs(yFlatBoundaryVector(end)-yFlatBoundaryVector(1)));

        reservedSize = 2 * (shiftHoloDist + zFlatBoundary(1)) * ...
            (0.5 * maxTransverseSize / (TransducerSf.radiusOfCurvature - zFlatBoundary(1)));
        
        xMin = min(xFlatBoundaryVector);
        xMax = max(xFlatBoundaryVector);
        yMin = min(yFlatBoundaryVector);
        yMax = max(yFlatBoundaryVector);

        xFlatBoundaryVectorResized = (xMin - reservedSize/2):dxFlatBoundary:(xMax + reservedSize/2);
        yFlatBoundaryVectorResized = (yMin - reservedSize/2):dyFlatBoundary:(yMax + reservedSize/2);

        % Use the same grid convention (ndgrid or meshgrid) as the input
        if isNdgrid
            [xFlatBoundary, yFlatBoundary] = ndgrid(xFlatBoundaryVectorResized, yFlatBoundaryVectorResized);
        else
            [xFlatBoundary, yFlatBoundary] = meshgrid(xFlatBoundaryVectorResized, yFlatBoundaryVectorResized);
        end
        zFlatBoundary = -shiftHoloDist * ones(size(xFlatBoundary));


        % Back‑project the source vibrational velocity to the shifted plane
        SourceParameters = [];
        SourceParameters.xGrid = xFlatBoundary;
        SourceParameters.yGrid = yFlatBoundary;
        SourceParameters.zGrid = zFlatBoundary;

        FieldParameters = [];
        FieldParameters.xGrid = TransducerSf.xGrid(activeTransducer);
        FieldParameters.yGrid = TransducerSf.yGrid(activeTransducer);
        FieldParameters.zGrid = TransducerSf.zGrid(activeTransducer);
        FieldParameters.dx = TransducerSf.dx;
        FieldParameters.dy = TransducerSf.dy;
        FieldParameters.input = TransducerSf.complexVelocityAmplitude(activeTransducer);

        regime = 6; % Back‑projection: V on a sphere --> P on planes
        isTransient = false;
        pFlatBoundary = rayleigh_simulator(expSign, frequency, regime, ...
            simulationDevice, isTransient, ...
            SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

        % Forward‑project from the shifted plane to the desired plane
        SourceParameters = [];
        SourceParameters.xGrid = xFlatBoundary;
        SourceParameters.yGrid = yFlatBoundary;
        SourceParameters.zGrid = zFlatBoundary;
        SourceParameters.dx = dxFlatBoundary;
        SourceParameters.dy = dyFlatBoundary;
        SourceParameters.input = pFlatBoundary;

        FieldParameters = [];
        FieldParameters.xGrid = xFlatBoundaryDesired;
        FieldParameters.yGrid = yFlatBoundaryDesired;
        FieldParameters.zGrid = zFlatBoundaryDesired;

        regime = 5; % Forward‑projection: P on a plane --> P at arbitrary points
        isTransient = false;
        pFlatBoundaryDesired = rayleigh_simulator(expSign, frequency, regime, ...
            simulationDevice, isTransient, ...
            SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

        % Return to the desired boundary condition at z = 0
        pFlatBoundary   = pFlatBoundaryDesired;

    else

        % Back‑project directly to the desired flat boundary condition
        SourceParameters = [];
        SourceParameters.xGrid = xFlatBoundary;
        SourceParameters.yGrid = yFlatBoundary;
        SourceParameters.zGrid = zFlatBoundary;

        FieldParameters = [];
        FieldParameters.xGrid = TransducerSf.xGrid(activeTransducer);
        FieldParameters.yGrid = TransducerSf.yGrid(activeTransducer);
        FieldParameters.zGrid = TransducerSf.zGrid(activeTransducer);
        FieldParameters.dx = TransducerSf.dx;
        FieldParameters.dy = TransducerSf.dy;
        FieldParameters.input = TransducerSf.complexVelocityAmplitude(activeTransducer);

        regime = 6; % Back‑projection: V on a sphere --> P on planes
        isTransient = false;
        pFlatBoundary = rayleigh_simulator(expSign, frequency, regime, ...
            simulationDevice, isTransient, ...
            SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);
    end

    % Convert pressure boundary condition into velocity on the flat surface
    vFlatBoundary = pressure_to_velocity_flat_surface(pFlatBoundary, ...
        dyFlatBoundary, dxFlatBoundary, frequency, Medium);
end


