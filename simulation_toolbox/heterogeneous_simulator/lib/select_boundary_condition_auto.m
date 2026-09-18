function izBoundaryCondition = select_boundary_condition_auto(MaterialMatrix, varargin)
%SELECT_BOUNDARY_CONDITION_AUTO  z-index for the boundary condition ('auto' mode).
%
%   Legacy (heterogeneous medium scan only):
%     iz = select_boundary_condition_auto(MaterialMatrix)
%     iz = select_boundary_condition_auto(MaterialMatrix, reserveIz)
%   Returns the last uniform xy-slice before the first non-uniform z-slice,
%   stepped back by reserveIz (default 5). If every slice is uniform, returns nzMedium.
%
%   Full routing (name-value pairs; use from xDDx_simulator-style scripts):
%     iz = select_boundary_condition_auto(MaterialMatrix, ...
%         'IsSpherical', isSphericalSource, ...
%         'WaterTest', waterTest, ...
%         'IzTarget', izTarget, ...
%         'RadiusOfCurvature', radiusOfCurvature, ...
%         'ReserveIz', 5)   % ReserveIz optional
%   Flat source  -> 1.
%   Spherical + water test (uniform medium) -> bowl apex index (same parity rule as izApex).
%   Spherical otherwise -> legacy uniform-slice scan above.

if isempty(varargin)
    izBoundaryCondition = iz_from_uniform_slice_scan(MaterialMatrix, 5);
    return
end

if numel(varargin) == 1 && isnumeric(varargin{1})
    izBoundaryCondition = iz_from_uniform_slice_scan(MaterialMatrix, varargin{1});
    return
end

p = inputParser;
p.addParameter('IsSpherical', [], @(x) isempty(x) || (isscalar(x) && (islogical(x) || isnumeric(x))));
p.addParameter('WaterTest', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
p.addParameter('IzTarget', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x == fix(x)));
p.addParameter('RadiusOfCurvature', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('ReserveIz', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.parse(varargin{:});
opts = p.Results;

if isempty(opts.IsSpherical)
    error('select_boundary_condition_auto: when using name-value arguments, IsSpherical must be specified.');
end
opts.IsSpherical = logical(opts.IsSpherical);

if ~opts.IsSpherical
    izBoundaryCondition = 1;
    return
end

if opts.WaterTest
    if isempty(opts.IzTarget) || isempty(opts.RadiusOfCurvature)
        error('select_boundary_condition_auto: WaterTest requires non-empty IzTarget and RadiusOfCurvature.');
    end
    dzMedium = MaterialMatrix.dz;
    izBoundaryCondition = opts.IzTarget - (opts.RadiusOfCurvature / dzMedium);
    if mod(izBoundaryCondition, 2) ~= 0
        izBoundaryCondition = fix(izBoundaryCondition) + 1;
    end
    nzMedium = size(MaterialMatrix.c0, 3);
    if izBoundaryCondition < 1
        izBoundaryCondition = 1;
    end
    if izBoundaryCondition > nzMedium
        error(['select_boundary_condition_auto: spherical water-test bowl apex z-index is %g, outside [1,%d].\n' ...
            'Increase izTarget or use a smaller radius of curvature.'], ...
            izBoundaryCondition, nzMedium);
    end
    return
end

izBoundaryCondition = iz_from_uniform_slice_scan(MaterialMatrix, opts.ReserveIz);
end


function izBoundaryCondition = iz_from_uniform_slice_scan(MaterialMatrix, reserveIz)
nzMedium = size(MaterialMatrix.c0, 3);
izBoundaryCondition = [];

for iz = 1 : nzMedium
    diffC0 = abs(MaterialMatrix.c0(:,:,iz) - MaterialMatrix.c0(1,1,iz));
    diffRho0 = abs(MaterialMatrix.rho0(:,:,iz) - MaterialMatrix.rho0(1,1,iz));
    if any(diffC0(:) > eps('single')) || any(diffRho0(:) > eps('single'))
        if iz == 1
            error('Auto boundary condition: the first z-slice (iz=1) is non-uniform; cannot select a uniform boundary position.');
        end
        izBoundaryCondition = iz - reserveIz;
        izBoundaryCondition = max(1, izBoundaryCondition);
        return;
    end
end

izBoundaryCondition = nzMedium;
end
