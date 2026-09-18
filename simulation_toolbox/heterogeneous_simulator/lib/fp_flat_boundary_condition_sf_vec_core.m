function [vFlatBoundary] = ...
    fp_flat_boundary_condition_sf_vec_core(TransducerSf, soundSpeed, density, simulationDevice, ...
    xFlatBoundary, yFlatBoundary, zFlatBoundary, dxFlatBoundary, dyFlatBoundary, ...
    ServiceParameters)

    if nargin < 10
        ServiceParameters = [];
        ServiceParameters.threadsPerBlockGPU = 128;
    end

    Medium = [];
    Medium.soundSpeed = soundSpeed;
    Medium.density = density;

% Ensure grids are 2D if a vector representation is provided
if isvector(TransducerSf.xGrid)
    TransducerSf = reshape_transducer_dim(TransducerSf);
end

TransducerSf = trim_TransducerSf_inactive_points(TransducerSf);

% Geometric and frequency parameters
expSign   = TransducerSf.expSign;
frequency = TransducerSf.frequency;

isSphericalSource = false;
if isfield(TransducerSf, 'radiusOfCurvature') && ~isempty(TransducerSf.radiusOfCurvature)
    radiusOfCurvature = TransducerSf.radiusOfCurvature;
    isSphericalSource = true;
else
    radiusOfCurvature = [];
end


% Forward‑project from the shifted plane to the desired plane z = 0
SourceParameters = [];
SourceParameters.xGrid = TransducerSf.xGrid;
SourceParameters.yGrid = TransducerSf.yGrid;
SourceParameters.zGrid = TransducerSf.zGrid;
SourceParameters.dx = TransducerSf.dx;
SourceParameters.dy = TransducerSf.dy;
SourceParameters.input = TransducerSf.complexVelocityAmplitude;

FieldParameters = [];
FieldParameters.xGrid = xFlatBoundary;
FieldParameters.yGrid = yFlatBoundary;
FieldParameters.zGrid = zFlatBoundary;

isTransient = false;

if isSphericalSource
    regime = 3; % 3 Forward-projection: V on a sphere --> P at an arbitrary set of points
    %Rayleigh simulator function with the complex pressure amplitude output
    [ pFlatBoundary ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

else

    regime = 4; % 4 Forward-projection: V on a plane --> P at an arbitrary set of points
    %Rayleigh simulator function with the complex pressure amplitude output
    [ pFlatBoundary ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters);

end

% Convert pressure boundary condition into velocity on the flat surface
vFlatBoundary = pressure_to_velocity_flat_surface(pFlatBoundary, dyFlatBoundary, dxFlatBoundary, frequency, Medium);

end
