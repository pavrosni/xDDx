function visualize_field(xSensor, ySensor, zSensor, complex_pressure, MaterialMatrix, ...
    xGridMedium, yGridMedium, zGridMedium, ixTarget, iyTarget, ...
    xBoundaryCondition, yBoundaryCondition, zBoundaryCondition, varargin)
%VISUALIZE_FIELD Visualize pressure field from k-Wave simulation.
%   visualize_field(xSensor, ySensor, zSensor, complex_pressure, MaterialMatrix, ...
%       xGridMedium, yGridMedium, zGridMedium, ixTarget, iyTarget)
%   visualize_field(..., Name, Value, ...)
%
%   This function visualizes the pressure field from a k-Wave simulation based on
%   the dimensionality of the data (3D, 2D, 1D, or 0D).
%
%   Inputs:
%       xSensor, ySensor, zSensor - Sensor coordinates (reshaped to 3D arrays)
%       complex_pressure          - Complex pressure field (reshaped to 3D array)
%       MaterialMatrix            - Material matrix structure with c0 field
%                                   or [] to skip medium overlay in 2D plots
%       xGridMedium, yGridMedium, zGridMedium - Medium grid coordinates
%       ixTarget, iyTarget        - Target indices for 2D slice extraction
%       xBoundaryCondition, yBoundaryCondition, zBoundaryCondition - Boundary condition coordinates
%
%   Optional Name-Value Pair Arguments:
%       'levelArray'              - Array of isosurface levels for 3D plot (default: [0.1 0.3 0.5 0.7 0.9])
%       'transparencyArray'       - Array of transparency values for each level (default: [0.1 0.2 0.3 0.4 0.5])
%       'pressureTransparency'    - Transparency for 2D overlay plots (default: 0.5)
%       'aperture'                - Transducer aperture in meters (default: 50e-3).
%                                   Pass [apertureX apertureY] for asymmetric sources.
%       'radiusOfCurvature'       - Transducer radius of curvature in meters (default: 50e-3)
%       'xDDxPath'                - Path to the xDDx toolbox library (default: '')
%
%   Example:
%       visualize_field(xSensor, ySensor, zSensor, complex_pressure, MaterialMatrix, ...
%           xGridMedium, yGridMedium, zGridMedium, ixTarget, iyTarget, ...
%           'levelArray', [0.1 0.3 0.5], 'pressureTransparency', 0.6);

% Parse optional arguments
p = inputParser;
p.addParameter('levelArray', [0.1 0.3 0.5 0.7 0.9], @isnumeric);
p.addParameter('transparencyArray', [0.1 0.2 0.3 0.4 0.5], @isnumeric);
p.addParameter('pressureTransparency', 0.5, @(x) isnumeric(x) && x >= 0 && x <= 1);
p.addParameter('aperture', 50e-3, @isnumeric);
p.addParameter('radiusOfCurvature', 50e-3, @isnumeric);
p.addParameter('xDDxPath', '', @ischar);
p.addParameter('transducerSf', [], @(x) isempty(x) || isstruct(x));

p.parse(varargin{:});
params = p.Results;

% Permute sensor data to create 3D field arrays
xField3D = permute(xSensor, [2 1 3]);
yField3D = permute(ySensor, [2 1 3]);
zField3D = permute(zSensor, [2 1 3]);
pField3D = permute(complex_pressure, [2 1 3]);

% Find maximum pressure and location
[pMax3D, iMax3D] = max(abs(pField3D(:)));
xMax3D = xField3D(iMax3D);
yMax3D = yField3D(iMax3D);
zMax3D = zField3D(iMax3D);

% Determine dimensionality and plot accordingly
if abs(ndims(squeeze(pField3D)) - 3) < eps
    % 3D field - plot isosurfaces
    plot_3d_field(xField3D, yField3D, zField3D, pField3D, pMax3D, ...
        xMax3D, yMax3D, zMax3D, params.levelArray, params.transparencyArray, ...
        params.xDDxPath, params.transducerSf);
    
elseif ~isvector(squeeze(pField3D))
    % 2D field - plot with or without medium overlay.
    if has_medium_overlay_data(MaterialMatrix, xGridMedium, yGridMedium, zGridMedium)
        plot_2d_field_with_overlay(xField3D, yField3D, zField3D, pField3D, ...
            MaterialMatrix, xGridMedium, yGridMedium, zGridMedium, ...
            ixTarget, iyTarget, ...
            xBoundaryCondition, yBoundaryCondition, zBoundaryCondition, ...
            params.pressureTransparency, ...
            params.aperture, params.radiusOfCurvature);
    else
        plot_2d_field(xField3D, yField3D, zField3D, pField3D);
    end
    
elseif numel(squeeze(pField3D)) > 1
    % 1D field
    [~, ~, ~, ~] = plot_1d_field(xField3D, yField3D, zField3D, pField3D);
    
else
    % 0D field (single point)
    [~, ~, ~, ~] = disp_0d_field(xField3D, yField3D, zField3D, pField3D);
end

end

function tf = has_medium_overlay_data(MaterialMatrix, xGridMedium, yGridMedium, zGridMedium)
tf = isstruct(MaterialMatrix) && isfield(MaterialMatrix, 'c0') && ~isempty(MaterialMatrix) && ...
    ~isempty(xGridMedium) && ~isempty(yGridMedium) && ~isempty(zGridMedium);
end
