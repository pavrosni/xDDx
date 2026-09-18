function plot_overlay_distributions(z, r, p, zCtLayer, rCtLayer, ctLayer, xCircle, yCircle, zBoundaryCondition, rBoundaryCondition, varargin)
%PLOT_OVERLAY_DISTRIBUTIONS Plot two distributions overlaid with different colormaps.
%   plot_overlay_distributions(z, r, p, zCtLayer, rCtLayer, ctLayer)
%   plot_overlay_distributions(z, r, p, zCtLayer, rCtLayer, ctLayer, Name, Value, ...)
%
%   This function plots two 2D distributions overlaid:
%   1. Background: ctLayer(zCtLayer, rCtLayer) with 'bone' colormap
%   2. Foreground: p(z, r) with 'jet' colormap (semi-transparent)
%
%   Inputs:
%       z, r, p          - Foreground distribution coordinates and values (plotted on top)
%       zCtLayer, rCtLayer, ctLayer - Background distribution coordinates and values (plotted below)
%       xCircle, yCircle - Circle coordinates (optional)
%
%   Optional Name-Value Pair Arguments:
%       'alpha'          - Transparency of foreground layer (0-1, default: 0.6)
%       'xlabel'         - X-axis label (default: 'z, mm')
%       'ylabel'         - Y-axis label (default: 'r, mm')
%       'title'          - Plot title (default: 'Overlaid Distributions')
%       'figHandle'      - Figure handle to plot in (default: new figure)
%       'axesHandle'     - Axes handle to plot in (default: new axes)
%       'backgroundCmap' - Colormap for background (default: 'bone')
%       'foregroundCmap' - Colormap for foreground (default: 'jet')
%       'colorbar'       - Show colorbar: 'both', 'background', 'foreground', 'none' (default: 'both')
%
%   Example:
%       plot_overlay_distributions(z, r, p, zCtLayer, rCtLayer, ctLayer, ...
%           'alpha', 0.7, 'title', 'Pressure over CT Layer');

% Parse input arguments
p_input = inputParser;
p_input.addParameter('alpha', 0.6, @(x) isnumeric(x) && x >= 0 && x <= 1);
p_input.addParameter('xlabel', 'z, mm', @ischar);
p_input.addParameter('ylabel', 'r, mm', @ischar);
p_input.addParameter('title', 'Overlaid Distributions', @ischar);
p_input.addParameter('figHandle', [], @(x) isempty(x) || isgraphics(x, 'figure'));
p_input.addParameter('axesHandle', [], @(x) isempty(x) || isgraphics(x, 'axes'));
p_input.addParameter('backgroundCmap', 'bone', @ischar);
p_input.addParameter('foregroundCmap', 'jet', @ischar);
p_input.addParameter('colorbar', 'none', @(x) ismember(x, {'both', 'background', 'foreground', 'none'}));

p_input.parse(varargin{:});
params = p_input.Results;

% Determine figure and axes
if isempty(params.figHandle)
    figHandle = figure;
else
    figHandle = params.figHandle;
    figure(figHandle);
end

if isempty(params.axesHandle)
    axesHandle = axes('Parent', figHandle);
else
    axesHandle = params.axesHandle;
    axes(axesHandle);
end

hold(axesHandle, 'on');

% Ensure inputs are 2D matrices
if ndims(p) > 2
    p = squeeze(p);
end
if ndims(ctLayer) > 2
    ctLayer = squeeze(ctLayer);
end
if ndims(z) > 2
    z = squeeze(z);
end
if ndims(r) > 2
    r = squeeze(r);
end
if ndims(zCtLayer) > 2
    zCtLayer = squeeze(zCtLayer);
end
if ndims(rCtLayer) > 2
    rCtLayer = squeeze(rCtLayer);
end

% Build coordinate grids for limit computation (significant area only)
if isvector(z) && isvector(r)
    [Z_fg_lim, R_fg_lim] = meshgrid(z, r);
else
    Z_fg_lim = z;
    R_fg_lim = r;
end
if isvector(zCtLayer) && isvector(rCtLayer)
    [Z_bg_lim, R_bg_lim] = meshgrid(zCtLayer, rCtLayer);
else
    Z_bg_lim = zCtLayer;
    R_bg_lim = rCtLayer;
end

% Determine axis limits from significant region only (where data exceeds 1% of max)
thresh_fg = 0.01 * max(p(:));
thresh_bg = 0.01 * max(ctLayer(:));
mask_fg = p >= thresh_fg;
mask_bg = ctLayer >= thresh_bg;
z_sig = [Z_fg_lim(mask_fg); Z_bg_lim(mask_bg)];
r_sig = [R_fg_lim(mask_fg); R_bg_lim(mask_bg)];

if ~isempty(z_sig)
    z_min = min(z_sig);
    z_max = max(z_sig);
    r_min = min(r_sig);
    r_max = max(r_sig);
    % Add small margin (2%) so plot is not cramped
    z_span = z_max - z_min;
    r_span = r_max - r_min;
    if z_span > 0
        z_min = z_min - 0.02 * z_span;
        z_max = z_max + 0.02 * z_span;
    end
    if r_span > 0
        r_min = r_min - 0.02 * r_span;
        r_max = r_max + 0.02 * r_span;
    end
else
    % Fallback: full data range
    z_min = min([min(z(:)), min(zCtLayer(:))]);
    z_max = max([max(z(:)), max(zCtLayer(:))]);
    r_min = min([min(r(:)), min(rCtLayer(:))]);
    r_max = max([max(r(:)), max(rCtLayer(:))]);
end

% Expand limits to include transducer arc (solid red) and boundary condition so they are always visible
if ~isempty(xCircle) && ~isempty(yCircle)
    z_min = min(z_min, min(xCircle(:)));
    z_max = max(z_max, max(xCircle(:)));
    r_min = min(r_min, min(yCircle(:)));
    r_max = max(r_max, max(yCircle(:)));
end
if ~isempty(zBoundaryCondition) && ~isempty(rBoundaryCondition)
    z_min = min(z_min, min(zBoundaryCondition(:)));
    z_max = max(z_max, max(zBoundaryCondition(:)));
    r_min = min(r_min, min(rBoundaryCondition(:)));
    r_max = max(r_max, max(rBoundaryCondition(:)));
end
% Re-apply margin after including overlay elements
z_span = z_max - z_min;
r_span = r_max - r_min;
if z_span > 0
    z_min = z_min - 0.02 * z_span;
    z_max = z_max + 0.02 * z_span;
end
if r_span > 0
    r_min = r_min - 0.02 * r_span;
    r_max = r_max + 0.02 * r_span;
end

% Plot background distribution (ctLayer)
% Use imagesc for better performance and automatic scaling.
% Handle both meshgrid-style and ndgrid-style coordinate matrices.
[z_bg_vec, r_bg_vec, ctLayerPlot] = get_plot_vectors_and_data(zCtLayer, rCtLayer, ctLayer);
h_bg = imagesc(axesHandle, z_bg_vec, r_bg_vec, ctLayerPlot);

set(h_bg, 'AlphaData', 1); % Fully opaque background
colormap(axesHandle, params.backgroundCmap);
caxis(axesHandle, caxis_increasing_limits(min(ctLayer(:)), max(ctLayer(:))));

% Set axis limits and properties for main axes
xlim(axesHandle, [z_min, z_max]);
ylim(axesHandle, [r_min, r_max]);
set(axesHandle, 'YDir', 'normal', 'Box', 'on');
axis(axesHandle, 'equal');

% Add colorbar for background if requested
if ismember(params.colorbar, {'both', 'background'})
    cb_bg = colorbar(axesHandle, 'southoutside');
    try
        cb_bg.Label.String = 'Medium Sound Speed, m/s';
    catch
        xlabel(cb_bg, 'Medium Sound Speed, m/s');
    end
end

% Plot foreground distribution (p) with transparency
% Create a new axes on top for the foreground with same position
% Set it up with no visible decorations from the start
axes2 = axes('Parent', figHandle, 'Position', get(axesHandle, 'Position'), ...
    'Color', 'none', 'Box', 'off', ...
    'XTick', [], 'YTick', [], ...
    'XColor', 'none', 'YColor', 'none');

if isvector(z) && isvector(r)
    z_fg_vec = z;
    r_fg_vec = r;
    pPlot = p;
else
    % If z and r are already 2D grids
    [z_fg_vec, r_fg_vec, pPlot] = get_plot_vectors_and_data(z, r, p);
end

h_fg = imagesc(axes2, z_fg_vec, r_fg_vec, pPlot);

hold(axes2, 'on');
% Plot circle if provided
if ~isempty(xCircle) && ~isempty(yCircle)
   plot(axes2, xCircle, yCircle, 'r-', 'LineWidth', 3);
end

if ~isempty(zBoundaryCondition) && ~isempty(rBoundaryCondition)
    orange = [0.8500 0.3250 0.0980];
    plot(axes2, zBoundaryCondition, rBoundaryCondition, '-', 'Color', orange, 'LineWidth', 2); 
end

function [xVec, yVec, dataOut] = get_plot_vectors_and_data(xGrid, yGrid, dataIn)
% Convert vectors or structured 2D grids to imagesc-compatible vectors/data.
if isvector(xGrid) && isvector(yGrid)
    xVec = xGrid(:).';
    yVec = yGrid(:);
    dataOut = dataIn;
    return;
end

% meshgrid-like: x varies along columns, y varies along rows
xDiffCols = abs(diff(xGrid, 1, 2));
xDiffRows = abs(diff(xGrid, 1, 1));
yDiffCols = abs(diff(yGrid, 1, 2));
yDiffRows = abs(diff(yGrid, 1, 1));
xVarCols = max(xDiffCols(:));
xVarRows = max(xDiffRows(:));
yVarCols = max(yDiffCols(:));
yVarRows = max(yDiffRows(:));

isMeshgridLike = (xVarCols >= xVarRows) && (yVarRows >= yVarCols);

if isMeshgridLike
    xVec = xGrid(1, :);
    yVec = yGrid(:, 1);
    dataOut = dataIn;
else
    % ndgrid-like: x varies along rows, y varies along columns
    xVec = xGrid(:, 1).';
    yVec = yGrid(1, :).';
    dataOut = dataIn.';
end
end

hold(axes2, 'off');



set(h_fg, 'AlphaData', params.alpha); % Semi-transparent foreground
colormap(axes2, params.foregroundCmap);
caxis(axes2, caxis_increasing_limits(min(p(:)), max(p(:))));

% Make axes2 background transparent and match limits
set(axes2, 'Color', 'none');
xlim(axes2, [z_min, z_max]);
ylim(axes2, [r_min, r_max]);
set(axes2, 'YDir', 'normal');
axis(axes2, 'equal');



% Set labels and title on the main axes
xlabel(axesHandle, params.xlabel);
ylabel(axesHandle, params.ylabel);
title(axesHandle, params.title);

% Link axes for zoom/pan (ensures both stay aligned)
linkaxes([axesHandle, axes2], 'xy');
% Link positions so axes2 stays over axesHandle on figure resize
hLinkPos = linkprop([axesHandle, axes2], 'Position');
setappdata(figHandle, 'PlotOverlayPositionLink', hLinkPos);

% Make sure axesHandle is the active axes (so it shows the box and labels)
%axes(axesHandle);

hold(axesHandle, 'off');

% Hide axes2 ticks, labels, and box completely (they'll show on axesHandle)
% Remove all axis decorations from axes2
set(axes2, 'XTick', [], 'YTick', [], ...
    'XTickLabel', [], 'YTickLabel', [], ...
    'XColor', 'none', 'YColor', 'none', ...
    'Box', 'off');
% Make sure the axes box is not drawn
set(axes2, 'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off');

% Add colorbar for foreground if requested
if ismember(params.colorbar, {'both', 'foreground'})
    cb_fg = colorbar(axes2);
    cb_fg.Label.String = 'Pressure amplitude, Pa';
    cb_fg.Location = 'eastoutside';
    
    % Re-hide axes decorations after colorbar (colorbar might reset some properties)
    set(axes2, 'XTick', [], 'YTick', [], ...
        'XTickLabel', [], 'YTickLabel', [], ...
        'XColor', 'none', 'YColor', 'none', ...
        'Box', 'off');
end

set(axes2,'visible','off');

% Re-apply limits to main axes (colorbar/axis equal may have reset axesHandle)
xlim(axesHandle, [z_min, z_max]);
ylim(axesHandle, [r_min, r_max]);

end


function lims = caxis_increasing_limits(lo, hi)
%CAXIS_INCREASING_LIMITS  [lo, hi] suitable for caxis (strictly increasing).
if hi > lo
    lims = [lo, hi];
    return
end
if ~isfinite(lo) || ~isfinite(hi)
    lims = [0, 1];
    return
end
scale = max(abs(lo), abs(hi));
delta = max(scale * eps, realmin('double'));
lims = [lo - delta, hi + delta];
end

