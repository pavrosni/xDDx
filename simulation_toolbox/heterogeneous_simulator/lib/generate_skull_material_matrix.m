function MaterialMatrix = generate_skull_material_matrix(varargin)
%GENERATE_SKULL_MATERIAL_MATRIX Generate MaterialMatrix for a human skull hemisphere model.
%   MaterialMatrix = generate_skull_material_matrix()
%   MaterialMatrix = generate_skull_material_matrix(Name, Value, ...)
%
%   This function generates a 3D MaterialMatrix structure representing a hemisphere
%   skull model with quasi-random bumps and dents on the inner surface. The skull
%   has uniform material properties throughout.
%
%   Optional Name-Value Pair Arguments:
%       'dx', 'dy', 'dz'          - Voxel size in meters (default: 0.5e-3 for all)
%       'c0'                      - Sound speed in m/s (default: 2331)
%       'rho0'                    - Density in kg/m^3 (default: 1732)
%       'alpha'                   - Attenuation in dB/cm (default: 8.83)
%       'outerRadius'             - Outer radius of hemisphere in meters (default: 0.09)
%       'skullThickness'          - Skull thickness in meters (default: 0.006)
%       'perturbationDepth'       - Maximum depth of bumps/dents in meters (default: 0.002)
%       'perturbationScale'       - Spatial scale of perturbations relative to radius (default: 0.15)
%       'gridSize'                - Grid size [Nx, Ny, Nz] (auto-calculated if not provided)
%       'water_c0'                - Sound speed of water/soft tissue in m/s (default: 1500)
%       'water_rho0'              - Density of water/soft tissue in kg/m^3 (default: 1000)
%       'water_alpha'             - Attenuation of water/soft tissue in dB/cm (default: 0)
%       'randomSeed'              - Random seed for reproducibility (default: random)
%
%   Output:
%       MaterialMatrix - Structure with fields:
%           .c0    - 3D matrix of sound speeds
%           .rho0  - 3D matrix of densities
%           .alpha - 3D matrix of attenuation coefficients
%           .dx, .dy, .dz - Grid spacing values
%
%   Example:
%       MaterialMatrix = generate_skull_material_matrix('outerRadius', 0.1, ...
%           'skullThickness', 0.007, 'perturbationDepth', 0.003);

% Parse input arguments
p = inputParser;
p.addParameter('dx', 0.5e-3, @(x) isnumeric(x) && x > 0);
p.addParameter('dy', 0.5e-3, @(x) isnumeric(x) && x > 0);
p.addParameter('dz', 0.5e-3, @(x) isnumeric(x) && x > 0);
p.addParameter('c0', 2331, @(x) isnumeric(x) && x > 0);
p.addParameter('rho0', 1732, @(x) isnumeric(x) && x > 0);
p.addParameter('alpha', 8.83, @(x) isnumeric(x) && x >= 0);
p.addParameter('outerRadius', 90e-3, @(x) isnumeric(x) && x > 0); % 90 mm, typical adult skull
p.addParameter('skullThickness', 6e-3, @(x) isnumeric(x) && x > 0); % 6 mm, typical skull thickness
p.addParameter('perturbationDepth', 4e-3, @(x) isnumeric(x) && x >= 0); % 3 mm max perturbation
p.addParameter('perturbationScale', 0.15, @(x) isnumeric(x) && x > 0); % Spatial scale of perturbations
p.addParameter('gridSize', [], @(x) isempty(x) || (isnumeric(x) && length(x) == 3 && all(x > 0)));
p.addParameter('water_c0', 1500, @(x) isnumeric(x) && x > 0);
p.addParameter('water_rho0', 1000, @(x) isnumeric(x) && x > 0);
p.addParameter('water_alpha', 0, @(x) isnumeric(x) && x >= 0);
p.addParameter('randomSeed', 42, @(x) isempty(x) || isnumeric(x));

p.parse(varargin{:});
params = p.Results;

% Set random seed if provided
if ~isempty(params.randomSeed)
    rng(params.randomSeed);
end

% Calculate grid size if not provided
if isempty(params.gridSize)
    % Add some padding around the hemisphere
    padding = params.outerRadius * 0.2; % 20% padding
    maxExtent = params.outerRadius + padding;
    Nx = ceil(2 * maxExtent / params.dx);
    Ny = ceil(2 * maxExtent / params.dy);
    Nz = ceil((params.outerRadius + padding) / params.dz);
    % Make dimensions odd for symmetry
    %Nx = Nx + mod(Nx+1, 2);
    %Ny = Ny + mod(Ny+1, 2);
    %Nz = Nz + mod(Nz+1, 2);
    Nx = 2*round(Nx/2);
    Ny = 2*round(Ny/2);
    Nz = 2*round(Nz/2);
else
    Nx = params.gridSize(1);
    Ny = params.gridSize(2);
    Nz = params.gridSize(3);
end

% Create coordinate grids
x = (-Nx/2:(Nx/2-1)) * params.dx;
y = (-Ny/2:(Ny/2-1)) * params.dy;
z = (0:(Nz-1)) * params.dz;
[X, Y, Z] = ndgrid(x, y, z);

% Calculate distance from origin in xy plane and total distance from origin
R_xy = sqrt(X.^2 + Y.^2);
R = sqrt(X.^2 + Y.^2 + Z.^2);

% Create smooth outer hemisphere surface (z >= 0)
% For a hemisphere: z = sqrt(outerRadius^2 - R_xy^2) for R_xy <= outerRadius
outerSurfaceZ = sqrt(params.outerRadius^2 - R_xy.^2);
outerSurfaceZ(R_xy > params.outerRadius) = 0;
outerSurfaceZ(Z < 0) = 0;

% Create inner surface with perturbations
% Use a combination of low-frequency noise for smooth bumps/dents
innerRadius = params.outerRadius - params.skullThickness;

% Generate quasi-random perturbations using Perlin-like noise
% Create multiple frequency components for natural-looking variations

scaleFreq = rand(3,1);
freqBase = 1 / (params.outerRadius * params.perturbationScale);

freqX = scaleFreq(1) * freqBase;
freqY = scaleFreq(2) * freqBase;
freqZ = scaleFreq(3) * freqBase;


amplitude = rand() * params.perturbationDepth;

    % Generate random phase shifts for each octave
phaseX = rand() * 2 * pi;
phaseY = rand() * 2 * pi;
phaseZ = rand() * 2 * pi;
    
% Create smooth perturbations using sinusoidal functions
perturbation = amplitude * sin(2*pi*freqX.*X + phaseX) .* ...
     sin(2*pi*freqY.*Y + phaseY) .* sin(2*pi*freqZ.*Z + phaseZ);

% Apply perturbations only to the inner surface region
% Normalize by distance from center to create radial variations
R_normalized = R_xy / params.outerRadius;
perturbation = perturbation .* (1 - R_normalized.^2); % Reduce perturbations near edges

% Create inner surface
innerSurfaceZ = sqrt(innerRadius^2 - R_xy.^2) + perturbation;
innerSurfaceZ(R_xy > innerRadius) = 0;
innerSurfaceZ(Z < 0) = 0;
% Ensure inner surface doesn't go below z=0 or above outer surface
innerSurfaceZ = max(0, min(innerSurfaceZ, outerSurfaceZ - params.skullThickness * 0.3));

% Create skull mask: points between inner and outer surfaces
% A point is in the skull if it's inside the outer hemisphere but outside the inner hemisphere
% Outer hemisphere: R <= outerRadius and Z >= 0 and Z <= outerSurfaceZ
% Inner hemisphere: Z >= 0 and Z <= innerSurfaceZ (accounting for perturbations)
insideOuterHemisphere = (R <= params.outerRadius) & (Z >= 0) & (Z <= outerSurfaceZ);
insideInnerHemisphere = (Z >= 0) & (Z <= innerSurfaceZ) & (R_xy <= innerRadius);
skullMask = insideOuterHemisphere & ~insideInnerHemisphere;

% Initialize material matrices with water properties
MaterialMatrix.c0 = params.water_c0 * ones(Nx, Ny, Nz);
MaterialMatrix.rho0 = params.water_rho0 * ones(Nx, Ny, Nz);
MaterialMatrix.alpha = params.water_alpha * ones(Nx, Ny, Nz);

% Assign skull properties to skull region
MaterialMatrix.c0(skullMask) = params.c0;
MaterialMatrix.rho0(skullMask) = params.rho0;
MaterialMatrix.alpha(skullMask) = params.alpha;

% Store grid spacing
MaterialMatrix.dx = params.dx;
MaterialMatrix.dy = params.dy;
MaterialMatrix.dz = params.dz;

MaterialMatrix.c0 = flip(MaterialMatrix.c0, 3);
MaterialMatrix.rho0 = flip(MaterialMatrix.rho0, 3);
MaterialMatrix.alpha = flip(MaterialMatrix.alpha, 3);

end