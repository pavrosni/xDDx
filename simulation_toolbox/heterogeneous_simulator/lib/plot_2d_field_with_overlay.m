function plot_2d_field_with_overlay(xField3D, yField3D, zField3D, pField3D, ...
    MaterialMatrix, xGridMedium, yGridMedium, zGridMedium, ...
    ixTarget, iyTarget, ...
    xBoundaryCondition, yBoundaryCondition, zBoundaryCondition, ...
    pressureTransparency, aperture, radiusOfCurvature)
%PLOT_2D_FIELD_WITH_OVERLAY Plot 2D field with CT layer overlay.

% Extract 2D field
[xField2D, yField2D, zField2D, pField2D] = plot_2d_field(xField3D, yField3D, zField3D, pField3D);
zBoundaryConditionRaw = zBoundaryCondition;
[xGridMediumVec, yGridMediumVec, zGridMediumVec, useVectorMediumGrid] = ...
    get_medium_grid_vectors(xGridMedium, yGridMedium, zGridMedium);

% Determine which plane (zy, xz, or xy) by checking constant coordinates.
isSphericalSource = ~isempty(radiusOfCurvature);
[apertureX, apertureY] = parse_aperture(aperture);
isXConstant = all(abs(xField2D(:) - xField2D(1)) < eps('single'));
isYConstant = all(abs(yField2D(:) - yField2D(1)) < eps('single'));
isZConstant = all(abs(zField2D(:) - zField2D(1)) < eps('single'));

if isXConstant && ~isZConstant
    % zy plane (x is constant)
    z = (squeeze(zField2D));
    r = (squeeze(yField2D));
    p = (squeeze(pField2D));
    
    ctLayer = squeeze(MaterialMatrix.c0(ixTarget,:,:));
    if useVectorMediumGrid
        zCtLayer = zGridMediumVec;
        rCtLayer = yGridMediumVec;
    else
        zCtLayer = squeeze(zGridMedium(1, :, :));
        rCtLayer = squeeze(yGridMedium(1, :, :));
    end
    
    if isSphericalSource
        zBoundaryCondition = [min(zBoundaryCondition(:)) max(zBoundaryCondition(:))];
        rBoundaryCondition = [min(yBoundaryCondition(:)) max(yBoundaryCondition(:))];
    else
        zBoundaryCondition = [];
        rBoundaryCondition = [];
    end

    rLabelPlot = 'y, mm';
    titlePlot = 'Pressure amplitude, Pa (zy plane)';
    apertureInPlane = apertureY;

elseif isYConstant && ~isZConstant
    % xz plane (y is constant)
    z = (squeeze(zField2D));
    r = (squeeze(xField2D));
    p = (squeeze(pField2D));
    
    ctLayer = squeeze(MaterialMatrix.c0(:,iyTarget,:));
    if useVectorMediumGrid
        zCtLayer = zGridMediumVec;
        rCtLayer = xGridMediumVec;
    else
        zCtLayer = squeeze(zGridMedium(:, 1, :));
        rCtLayer = squeeze(xGridMedium(:, 1, :));
    end
    
    if isSphericalSource
        zBoundaryCondition = [min(zBoundaryCondition(:)) max(zBoundaryCondition(:))];
        rBoundaryCondition = [min(xBoundaryCondition(:)) max(xBoundaryCondition(:))];
    else
        zBoundaryCondition = [];
        rBoundaryCondition = [];
    end

    rLabelPlot = 'x, mm';
    titlePlot = 'Pressure amplitude, Pa (xz plane)';
    apertureInPlane = apertureX;

elseif isZConstant
    % xy plane (z is constant)
    z = squeeze(xField2D);
    r = squeeze(yField2D);
    p = squeeze(pField2D);

    zSliceValue = zField2D(1);
    if useVectorMediumGrid
        [~, izTarget] = min(abs(zGridMediumVec - zSliceValue));
    else
        [~, izTarget] = min(abs(squeeze(zGridMedium(1,1,:)) - zSliceValue));
    end
    ctLayer = squeeze(MaterialMatrix.c0(:,:,izTarget));
    if useVectorMediumGrid
        zCtLayer = xGridMediumVec;
        rCtLayer = yGridMediumVec;
    else
        zCtLayer = squeeze(xGridMedium(:, :, izTarget));
        rCtLayer = squeeze(yGridMedium(:, :, izTarget));
    end

    % Do not draw transducer and boundary overlays for xy slices.
    zBoundaryCondition = [];
    rBoundaryCondition = [];
    xCircle = [];
    yCircle = [];

    rLabelPlot = 'y, mm';
    titlePlot = 'Pressure amplitude, Pa (xy plane)';
    apertureInPlane = max(apertureX, apertureY);

else
    error('Unable to determine 2D plane orientation from field grids.');
end

% Calculate transducer overlay.
if ~isZConstant && isSphericalSource
    angle1 = pi - asin(apertureInPlane/2/radiusOfCurvature);
    angle2 = pi + asin(apertureInPlane/2/radiusOfCurvature);
    circleParam = linspace(angle1, angle2, 100);
    xCircle = radiusOfCurvature*cos(circleParam) + radiusOfCurvature;
    yCircle = radiusOfCurvature*sin(circleParam);
elseif ~isZConstant
    % Flat transducer: draw as a solid line segment at the transducer plane.
    zFinite = zBoundaryConditionRaw(isfinite(zBoundaryConditionRaw));
    if isempty(zFinite)
        zTransducer = min(z(:));
    else
        zTransducer = mean(zFinite(:));
    end
    xCircle = [zTransducer, zTransducer];
    yCircle = [-apertureInPlane/2, apertureInPlane/2];
end

% Plot overlay distributions in millimeters to match the axis labels.
meterToMillimeter = 1e3;
plot_overlay_distributions(meterToMillimeter*z, meterToMillimeter*r, abs(p), ...
    meterToMillimeter*zCtLayer, meterToMillimeter*rCtLayer, ctLayer, ...
    meterToMillimeter*xCircle, meterToMillimeter*yCircle, ...
    meterToMillimeter*zBoundaryCondition, meterToMillimeter*rBoundaryCondition,...
    'alpha', pressureTransparency, ...
    'title', titlePlot, ...
    'xlabel', ternary_label(isZConstant, 'x, mm', 'z, mm'), ...
    'ylabel', rLabelPlot, ...
    'colorbar', 'both');

end

function out = ternary_label(condition, trueValue, falseValue)
if condition
    out = trueValue;
else
    out = falseValue;
end
end

function [xVec, yVec, zVec, useVectorGrid] = get_medium_grid_vectors(xGridMedium, yGridMedium, zGridMedium)
useVectorGrid = isvector(xGridMedium) && isvector(yGridMedium) && isvector(zGridMedium);

if useVectorGrid
    xVec = xGridMedium(:);
    yVec = yGridMedium(:);
    zVec = zGridMedium(:);
else
    xVec = [];
    yVec = [];
    zVec = [];
end
end

function [apertureX, apertureY] = parse_aperture(aperture)
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
end
