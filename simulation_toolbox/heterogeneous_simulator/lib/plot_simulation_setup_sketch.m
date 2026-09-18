function plot_simulation_setup_sketch(ax, plane, MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, aperture, radiusOfCurvature, xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, xGridSimulationVec, yGridSimulationVec, zGridSimulationVec, showSimBox, radialReserveX, radialReserveY, showLegend, centerLegend)
%PLOT_SIMULATION_SETUP_SKETCH Plot setup sketch in a given axis.
% plane: 'yz' (x fixed through target) or 'xz' (y fixed through target)
% showSimBox (optional): if false, the simulation box (blue) is not drawn; default true.
% radialReserveX, radialReserveY (optional): when provided, the boundary
%   condition transverse extents match generate_simulation_model; otherwise full uSim extents are used.
% showLegend (optional): if false, no legend is drawn; default true.
% centerLegend (optional): if true, center the legend under the figure; default false.

if nargin < 19
    showSimBox = true;
end
if nargin < 20
    radialReserveX = [];
end
if nargin < 21
    radialReserveY = [];
end
if nargin < 22
    showLegend = true;
end
if nargin < 23
    centerLegend = false;
end

if ~(isgraphics(ax, 'axes') && isscalar(ax))
    error('First argument must be a scalar axes handle.');
end
if isstring(plane)
    plane = char(plane);
end
plane = lower(plane);

dx = MaterialMatrix.dx;
dy = MaterialMatrix.dy;
dz = MaterialMatrix.dz;
isSphericalSource = ~isempty(radiusOfCurvature);
[apertureX, apertureY, apertureScalar] = parse_aperture(aperture);

nx = size(MaterialMatrix.c0, 1);
ny = size(MaterialMatrix.c0, 2);
nz = size(MaterialMatrix.c0, 3);

xVec = ((1:nx) - ixTarget) .* dx;              % meters, target at 0
yVec = ((1:ny) - iyTarget) .* dy;              % meters, target at 0
if isSphericalSource
    zVec = ((1:nz) - izTarget) .* dz + radiusOfCurvature; % meters, target at z=R
    zBC  = ((izBoundaryCondition - izTarget) .* dz + radiusOfCurvature); % meters
    zTargetCoord = radiusOfCurvature;
else
    zVec = ((1:nz) - izBoundaryCondition) .* dz; % meters, boundary at z=0
    zBC  = 0;
    zTargetCoord = (izTarget - izBoundaryCondition) .* dz;
end

if nargin < 17 || isempty(xGridSimulationVec) || isempty(yGridSimulationVec) || isempty(zGridSimulationVec)
    xGridSimulationVec = [min(xVec) max(xVec)];
    yGridSimulationVec = [min(yVec) max(yVec)];
    zGridSimulationVec = [min(zVec) max(zVec)];
end

uSimX = [min(xGridSimulationVec(:)) max(xGridSimulationVec(:))]; % meters
uSimY = [min(yGridSimulationVec(:)) max(yGridSimulationVec(:))]; % meters
zSim  = [min(zGridSimulationVec(:)) max(zGridSimulationVec(:))]; % meters

% Boundary transverse extent: same formula as generate_simulation_model when reserves are given
if ~isempty(radialReserveX) && ~isempty(radialReserveY)
    if isSphericalSource
    zBCm = (izBoundaryCondition - izTarget) * dz + radiusOfCurvature;  % z at boundary in m (same as zGridMediumVec(izBC))
    transverseScale = (radiusOfCurvature - zBCm) / sqrt(radiusOfCurvature^2 - (apertureScalar/2)^2);
    transverseSizeBoundaryConditionX = (1 + radialReserveX) * apertureX * transverseScale;
    transverseSizeBoundaryConditionY = (1 + radialReserveY) * apertureY * transverseScale;
    else
    transverseSizeBoundaryConditionX = (1 + radialReserveX) * apertureX;
    transverseSizeBoundaryConditionY = (1 + radialReserveY) * apertureY;
    end
    transverseSizeBoundaryConditionX = max(transverseSizeBoundaryConditionX, 2*dx);  % at least 2 grid steps
    transverseSizeBoundaryConditionY = max(transverseSizeBoundaryConditionY, 2*dy);  % at least 2 grid steps
    xBoundaryHalf = transverseSizeBoundaryConditionX / 2;  % meters
    yBoundaryHalf = transverseSizeBoundaryConditionY / 2;  % meters
else
    xBoundaryHalf = [];  % use full uSim extent
    yBoundaryHalf = [];  % use full uSim extent
end

hold(ax, 'on');
set(ax, 'YDir', 'normal', 'Box', 'on');
axis(ax, 'equal');

flatTransducerZmm = [];
flatTransducerUlimsMm = [];
switch plane
    case 'yz'
        bg = squeeze(MaterialMatrix.c0(ixTarget, :, :)); % [ny x nz]
        imagesc(ax, zVec * 1e3, yVec * 1e3, bg);
        colormap(ax, bone);
        xlabel(ax, '\it z\rm, mm', 'Interpreter', 'tex');
        ylabel(ax, '\it y\rm, mm', 'Interpreter', 'tex');
        title(ax, '\it yz\rm plane through target', 'Interpreter', 'tex');
        if ~isempty(yBoundaryHalf)
            yBoundaryLims = [-yBoundaryHalf, yBoundaryHalf] * 1e3;  % mm
        else
            yBoundaryLims = uSimY * 1e3;
        end
        yTransducerLims = [-apertureY/2, apertureY/2] * 1e3;  % mm
        [yBoxMin, yBoxMax] = clip_interval(yFieldBegin, yFieldEnd, uSimY(1), uSimY(2));
        [zBoxMin, zBoxMax] = clip_interval(zFieldBegin, zFieldEnd, zSim(1), zSim(2));
        boxMin = [zBoxMin, yBoxMin] * 1e3;
        boxSize = [zBoxMax - zBoxMin, yBoxMax - yBoxMin] * 1e3;
        outputFieldIsPoint = (boxSize(1) <= 1e-12 && boxSize(2) <= 1e-12);
        zBoxMinMm = boxMin(1);
        zBoxMaxMm = boxMin(1) + boxSize(1);
        targetPt = [zTargetCoord, 0] * 1e3;          % [z, y] in mm
        if isSphericalSource
            rim1 = [radiusOfCurvature - radiusOfCurvature*cos(asin((apertureY/2)/radiusOfCurvature)), +(apertureY/2)] * 1e3;
            rim2 = [radiusOfCurvature - radiusOfCurvature*cos(asin((apertureY/2)/radiusOfCurvature)), -(apertureY/2)] * 1e3;
            plotTransducerArc(ax, apertureY, radiusOfCurvature);
            plotFocusingRays(ax, targetPt, rim1, rim2);
        else
            plotTransducerLine(ax, zBC * 1e3, yTransducerLims);
            flatTransducerZmm = zBC * 1e3;
            flatTransducerUlimsMm = yTransducerLims;
        end
        plotTarget(ax, targetPt);
        if showSimBox
            plotSimBox(ax, [zSim(1), uSimY(1)] * 1e3, [(zSim(2)-zSim(1)), (uSimY(2)-uSimY(1))] * 1e3);  % Simulation box (full grid)
        end
        plotOutputFieldBox(ax, boxMin, boxSize);  % Output field region
        if isSphericalSource
            hBoundary = plotBoundary(ax, zBC * 1e3, yBoundaryLims);
            set(hBoundary, 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none');
        end
        setViewLimits(ax, yVec([1 end]), uSimY, zVec([1 end]), zSim, zBC);

    case 'xz'
        bg = squeeze(MaterialMatrix.c0(:, iyTarget, :)); % [nx x nz]
        imagesc(ax, zVec * 1e3, xVec * 1e3, bg);
        colormap(ax, bone);
        xlabel(ax, '\it z\rm, mm', 'Interpreter', 'tex');
        ylabel(ax, '\it x\rm, mm', 'Interpreter', 'tex');
        title(ax, '\it xz\rm plane through target', 'Interpreter', 'tex');
        if ~isempty(xBoundaryHalf)
            xBoundaryLims = [-xBoundaryHalf, xBoundaryHalf] * 1e3;  % mm
        else
            xBoundaryLims = uSimX * 1e3;
        end
        xTransducerLims = [-apertureX/2, apertureX/2] * 1e3;  % mm
        [xBoxMin, xBoxMax] = clip_interval(xFieldBegin, xFieldEnd, uSimX(1), uSimX(2));
        [zBoxMin, zBoxMax] = clip_interval(zFieldBegin, zFieldEnd, zSim(1), zSim(2));
        boxMin = [zBoxMin, xBoxMin] * 1e3;
        boxSize = [zBoxMax - zBoxMin, xBoxMax - xBoxMin] * 1e3;
        outputFieldIsPoint = (boxSize(1) <= 1e-12 && boxSize(2) <= 1e-12);
        zBoxMinMm = boxMin(1);
        zBoxMaxMm = boxMin(1) + boxSize(1);
        targetPt = [zTargetCoord, 0] * 1e3;          % [z, x] in mm
        if isSphericalSource
            rim1 = [radiusOfCurvature - radiusOfCurvature*cos(asin((apertureX/2)/radiusOfCurvature)), +(apertureX/2)] * 1e3;
            rim2 = [radiusOfCurvature - radiusOfCurvature*cos(asin((apertureX/2)/radiusOfCurvature)), -(apertureX/2)] * 1e3;
            plotTransducerArc(ax, apertureX, radiusOfCurvature);
            plotFocusingRays(ax, targetPt, rim1, rim2);
        else
            plotTransducerLine(ax, zBC * 1e3, xTransducerLims);
            flatTransducerZmm = zBC * 1e3;
            flatTransducerUlimsMm = xTransducerLims;
        end
        plotTarget(ax, targetPt);
        if showSimBox
            plotSimBox(ax, [zSim(1), uSimX(1)] * 1e3, [(zSim(2)-zSim(1)), (uSimX(2)-uSimX(1))] * 1e3);  % Simulation box (full grid)
        end
        plotOutputFieldBox(ax, boxMin, boxSize);  % Output field region
        if isSphericalSource
            hBoundary = plotBoundary(ax, zBC * 1e3, xBoundaryLims);
            set(hBoundary, 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none');
        end
        setViewLimits(ax, xVec([1 end]), uSimX, zVec([1 end]), zSim, zBC);

    otherwise
        error('Unknown plane "%s". Expected ''yz'' or ''xz''.', plane);
end

% Force flat-transducer line to be drawn last so it stays in front of the
% simulation box and output-region overlays.
if ~isSphericalSource && ~isempty(flatTransducerUlimsMm)
    plotTransducerLine(ax, flatTransducerZmm, flatTransducerUlimsMm);
end

axPosFixed = get(ax, 'Position');

% Use proxy handles for compatibility with older MATLAB legend limitations.
hTransLegend = plot(ax, nan, nan, 'r-', 'LineWidth', 2);
if isSphericalSource
    hBoundaryLegend = plot(ax, nan, nan, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
end
if showSimBox
    hSimBoxLegend = plot(ax, nan, nan, '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
end
if outputFieldIsPoint
    hFieldBoxLegend = plot(ax, nan, nan, 's', 'Color', [0 0.6 0.2], 'MarkerSize', 8, 'LineWidth', 1.5);
else
    hFieldBoxLegend = plot(ax, nan, nan, '--', 'Color', [0 0.6 0.2], 'LineWidth', 1.5);
end

cb = colorbar(ax, 'southoutside');
try
    cb.Label.String = 'Medium Sound Speed, m/s';
catch
    xlabel(cb, 'Medium Sound Speed, m/s');
end

lgd = [];
if showLegend
    if showSimBox
        if isSphericalSource
            lgd = legend(ax, [hTransLegend, hBoundaryLegend, hSimBoxLegend, hFieldBoxLegend], ...
                {'Transducer', 'Boundary condition', 'Simulation box', 'Output field region'});
        else
            lgd = legend(ax, [hTransLegend, hSimBoxLegend, hFieldBoxLegend], ...
                {'Transducer', 'Simulation box', 'Output field region'});
        end
    else
        if isSphericalSource
            lgd = legend(ax, [hTransLegend, hBoundaryLegend, hFieldBoxLegend], ...
                {'Transducer', 'Boundary condition', 'Output field region'});
        else
            lgd = legend(ax, [hTransLegend, hFieldBoxLegend], ...
                {'Transducer', 'Output field region'});
        end
    end
    lgd.Location = 'southoutside';
    lgd.Box = 'off';
    try
        lgd.Orientation = 'horizontal';
    catch
    end
    if centerLegend
        try
            lgd.NumColumns = 2;
        catch
        end
    end
end

% Restore axes position after MATLAB southoutside auto-layout shifts it.
set(ax, 'Position', axPosFixed);
stack_colorbar_and_legend(ax, axPosFixed, cb, lgd, zBoxMinMm, zBoxMaxMm, centerLegend);

hold(ax, 'off');

end

function h = plotTransducerArc(ax, aperture, R)
thetaMax = asin((aperture/2) / R);
theta = linspace(-thetaMax, thetaMax, 300);
z = R - R * cos(theta);
u = R * sin(theta); % y (for yz) or x (for xz)
h = plot(ax, z * 1e3, u * 1e3, 'r-', 'LineWidth', 2);
end

function h = plotTransducerLine(ax, zPosMm, uLimsMm)
h = line(ax, [zPosMm zPosMm], uLimsMm, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none');
end

function plotFocusingRays(ax, targetPt, rim1, rim2)
plot(ax, [rim1(1), targetPt(1)], [rim1(2), targetPt(2)], 'r--', 'LineWidth', 1.5);
plot(ax, [rim2(1), targetPt(1)], [rim2(2), targetPt(2)], 'r--', 'LineWidth', 1.5);
end

function plotTarget(ax, targetPt)
plot(ax, targetPt(1), targetPt(2), 'ro', 'LineWidth', 1.8, 'MarkerSize', 7);
end

function h = plotBoundary(ax, zBCmm, uLimsMm)
orange = [0.8500 0.3250 0.0980];
h = line(ax, [zBCmm zBCmm], uLimsMm, ...
    'Color', orange, 'LineStyle', '-', 'LineWidth', 2, 'Marker', 'none');
end

function h = plotSimBox(ax, boxMin, boxSize)
blue = [0 0.4470 0.7410];
h = rectangle(ax, 'Position', [boxMin(1) boxMin(2) boxSize(1) boxSize(2)], ...
    'EdgeColor', blue, 'LineStyle', '--', 'LineWidth', 2);
end

function h = plotOutputFieldBox(ax, boxMin, boxSize)
green = [0 0.6 0.2];
% Output field region can be a single point (xFieldBegin==xFieldEnd, etc.)
tol = 1e-12;
if boxSize(1) <= tol && boxSize(2) <= tol
    % Single point: plot a square marker at the output field location
    h = plot(ax, boxMin(1), boxMin(2), 's', 'Color', green, 'MarkerSize', 8, 'LineWidth', 1.5);
else
    h = rectangle(ax, 'Position', [boxMin(1) boxMin(2) boxSize(1) boxSize(2)], ...
        'EdgeColor', green, 'LineStyle', '--', 'LineWidth', 1.5);
end
end

function setViewLimits(ax, uMediumEnds, uSimEnds, zMediumEnds, zSimEnds, zBC)
uMin0 = min([uMediumEnds(:); uSimEnds(:)]);
uMax0 = max([uMediumEnds(:); uSimEnds(:)]);
uPad = 0.05 * (uMax0 - uMin0);
if uPad <= 0
    uPad = 0.01 * max([abs(uMax0), abs(uMin0), eps]);
end

% Always include z=0 where the transducer arc is defined.
zMin0 = min([zMediumEnds(:); zSimEnds(:); zBC; 0]);
zMax0 = max([zMediumEnds(:); zSimEnds(:); zBC]);
zPad = 0.05 * (zMax0 - zMin0);
if zPad <= 0
    zPad = 0.01 * max([abs(zMax0), abs(zMin0), eps]);
end

xlim(ax, [(zMin0 - zPad) (zMax0 + zPad)] * 1e3);
ylim(ax, [(uMin0 - uPad) (uMax0 + uPad)] * 1e3);
end

function [vMinOut, vMaxOut] = clip_interval(vMinIn, vMaxIn, clipMin, clipMax)
vMin = min(vMinIn, vMaxIn);
vMax = max(vMinIn, vMaxIn);
vMinOut = max(vMin, clipMin);
vMaxOut = min(vMax, clipMax);
if vMinOut > vMaxOut
    vMinOut = clipMin;
    vMaxOut = clipMin;
end
end

function stack_colorbar_and_legend(ax, axPos, cb, lgd, zBoxMinMm, zBoxMaxMm, centerLegend)
% Place colorbar below xlabel and legend below colorbar.
xL = xlim(ax);
xSpan = xL(2) - xL(1);
if xSpan <= 0
    fracMin = 0;
    fracMax = 1;
else
    fracMin = (zBoxMinMm - xL(1)) / xSpan;
    fracMax = (zBoxMaxMm - xL(1)) / xSpan;
end
fracMin = max(0, min(1, fracMin));
fracMax = max(0, min(1, fracMax));
if fracMax <= fracMin
    fracMin = 0;
    fracMax = 1;
end

cbW = axPos(3) * (fracMax - fracMin);
% Keep width tied to simulation-box span but center the stack under axes.
cbX = axPos(1) + 0.5 * (axPos(3) - cbW);
cbH = 0.020;
lgdH = 0.028;
gapCB = 0.065;
gapLG = 0.075;

% Keep colorbar/legend clearly below x-axis label and above window controls.
cbY = axPos(2) - cbH - gapCB;
lgdY = cbY - lgdH - gapLG;

set(cb, 'Units', 'normalized');
set(cb, 'Position', [cbX, cbY, cbW, cbH]);

if isempty(lgd) || ~isgraphics(lgd)
    return;
end

set(lgd, 'Units', 'normalized');
drawnow;

% Use the legend's natural width so the entries are not left-aligned inside
% an oversized legend rectangle.
lgdPos = get(lgd, 'Position');
lgdW = lgdPos(3);
if centerLegend
    lgdX = 0.5 - 0.5 * lgdW;
else
    lgdW = min(axPos(3), max(0.60 * axPos(3), cbW + 0.08));
    lgdX = axPos(1) + 0.5 * (axPos(3) - lgdW);
end
set(lgd, 'Position', [lgdX, lgdY, lgdW, lgdH]);
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
