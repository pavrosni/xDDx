function MaterialMatrix = load_material_matrix(dataPath)
%LOAD_MATERIAL_MATRIX Load and validate MaterialMatrix from a MAT file.
%
%   MaterialMatrix = load_material_matrix(dataPath)
%
% INPUT DATA FORMAT WITH HETEROGENEOUS MEDIUM PARAMETERS: a MAT file
% 'dataPath' with the following variable
%
% 'MaterialMatrix' scalar struct with fields:
%   'c0': 3D matrix with the sound speed in m/s at each grid node
%       of the heterogeneous medium (Cartesian grid only)
%   'rho0': 3D matrix with the density in kg/m^3 at each grid node
%       of the heterogeneous medium
%   'alpha': 3D matrix with the attenuation in dB/cm at the operating
%       frequency of the transducer (TransducerSf.frequency, given in Hz)
%       at each grid node of the heterogeneous medium. The simulator
%       converts these values to power law coefficients by dividing by
%       (TransducerSf.frequency * 1e-6)^alphaPower, where alphaPower is
%       specified separately in the simulator inputs
%   'dx': x-step of the medium Cartesian grid in m (positive scalar)
%   'dy': y-step of the medium Cartesian grid in m (positive scalar)
%   'dz': z-step of the medium Cartesian grid in m (positive scalar)
%
% All fields are required and must be non-empty and numeric. The matrices
% 'c0', 'rho0', and 'alpha' must have identical sizes [Nx, Ny, Nz], with
% the first, second, and third dimensions corresponding to x, y, and z,
% respectively (ndgrid format).
%
% Here, x and y are the transverse coordinates, and z is the coordinate
% along the direction of the beam.

    if nargin < 1 || isempty(dataPath)
        error('dataPath must be a non-empty char vector or string.');
    end
    if isa(dataPath, 'string')
        dataPath = char(dataPath);
    end
    if ~ischar(dataPath)
        error('dataPath must be a char vector or string.');
    end

    if exist(dataPath, 'file') ~= 2
        error('File not found: %s', dataPath);
    end

    S = load(dataPath);
    if ~isfield(S, 'MaterialMatrix')
        vars = strjoin(fieldnames(S), ', ');
        error('MAT file must contain variable ''MaterialMatrix''. Found: %s', vars);
    end

    MaterialMatrix = validate_material_matrix_structure(S.MaterialMatrix);
end

function MaterialMatrix = validate_material_matrix_structure(M)
    if ~isstruct(M) || numel(M) ~= 1
        error('MaterialMatrix must be a scalar struct.');
    end

    required = {'c0', 'rho0', 'alpha', 'dx', 'dy', 'dz'};
    for k = 1:numel(required)
        if ~isfield(M, required{k}) || isempty(M.(required{k}))
            error('MaterialMatrix missing or empty required field: %s', required{k});
        end
    end

    c0 = M.c0;
    rho0 = M.rho0;
    alpha = M.alpha;
    if ~isnumeric(c0) || ~isnumeric(rho0) || ~isnumeric(alpha)
        error('MaterialMatrix.c0, rho0, and alpha must be numeric arrays.');
    end
    if ndims(c0) ~= 3 || ndims(rho0) ~= 3 || ndims(alpha) ~= 3
        error('MaterialMatrix.c0, rho0, and alpha must be 3D arrays.');
    end
    if ~isequal(size(c0), size(rho0), size(alpha))
        error('MaterialMatrix.c0, rho0, and alpha must have identical sizes.');
    end

    dx = M.dx;
    dy = M.dy;
    dz = M.dz;
    if ~isnumeric(dx) || ~isscalar(dx) || dx <= 0
        error('MaterialMatrix.dx must be a positive scalar.');
    end
    if ~isnumeric(dy) || ~isscalar(dy) || dy <= 0
        error('MaterialMatrix.dy must be a positive scalar.');
    end
    if ~isnumeric(dz) || ~isscalar(dz) || dz <= 0
        error('MaterialMatrix.dz must be a positive scalar.');
    end

    MaterialMatrix = [];
    MaterialMatrix.c0 = c0;
    MaterialMatrix.rho0 = rho0;
    MaterialMatrix.alpha = alpha;
    MaterialMatrix.dx = double(dx);
    MaterialMatrix.dy = double(dy);
    MaterialMatrix.dz = double(dz);
end
