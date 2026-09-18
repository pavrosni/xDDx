function plot_pressure_with_ctlayer(z, r, p, ctLayer, zCtLayer, rCtLayer, varargin)
%PLOT_PRESSURE_WITH_CTLAYER Plot pressure distribution overlaid on ctLayer background
%
%   plot_pressure_with_ctlayer(z, y, p, ctLayer, zCtLayer, yCtLayer)
%   plot_pressure_with_ctlayer(..., 'Title', titleStr)
%   plot_pressure_with_ctlayer(..., 'XLabel', xLabelStr)
%   plot_pressure_with_ctlayer(..., 'YLabel', yLabelStr)
%
%   This function plots the pressure distribution (abs(p)) overlaid on a
%   background ctLayer distribution. The ctLayer is displayed with 'bone'
%   colormap, and the pressure is displayed with 'jet' colormap. The
%   pressure distribution prevails with variable transparency (higher
%   pressure = more opaque).
%
%   Inputs:
%       z - z coordinates for pressure (vector or matrix)
%       y - y coordinates for pressure (vector or matrix)
%       p - pressure data (matrix, same size as z and y grids)
%       ctLayer - ctLayer data matrix [Ny x (Nz-1)]
%       zCtLayer - z coordinates for ctLayer (vector)
%       yCtLayer - y coordinates for ctLayer (vector)
%       Optional name-value pairs:
%           'Title' - plot title (default: 'Pressure amplitude, Pa')
%           'XLabel' - x-axis label (default: 'z, mm')
%           'YLabel' - y-axis label (default: 'y, mm')
%
%   Example:
%       plot_pressure_with_ctlayer(z, y, p, ctLayer, zCtLayer, yCtLayer);

    % Parse optional inputs
    parser = inputParser;
    addParameter(parser, 'Title', 'Pressure amplitude, Pa', @ischar);
    addParameter(parser, 'XLabel', 'z, mm', @ischar);
    addParameter(parser, 'YLabel', 'y, mm', @ischar);
    addParameter(parser, 'PressureTransparency', 0.5, @isnumeric);
    parse(parser, varargin{:});
    
    titleStr = parser.Results.Title;
    xLabelStr = parser.Results.XLabel;
    yLabelStr = parser.Results.YLabel;
    pressureTransparency = parser.Results.PressureTransparency;
   
    ctLayer = interp2(zCtLayer, rCtLayer, ctLayer, z, r, 'linear', NaN);
    
    
    % Convert ctLayer to RGB using 'bone' colormap
    ctValid = ctLayer(~isnan(ctLayer(:)));
    if ~isempty(ctValid)
        ctMin = min(ctValid(:));
        ctMax = max(ctValid(:));
        if ctMax > ctMin
            ctNormalized = (ctLayer - ctMin) / (ctMax - ctMin);
        else
            ctNormalized = zeros(size(ctLayer));
        end
    else
        ctNormalized = zeros(size(ctLayer));
    end
    ctNormalized(isnan(ctNormalized)) = 0;
    boneMap = colormap('bone');
    ctRGB = ind2rgb(round(ctNormalized * (size(boneMap, 1) - 1)) + 1, boneMap);
    
    % Convert pressure to RGB using 'jet' colormap
    pAbs = abs(p);
    pValid = pAbs(~isnan(pAbs(:)));
    if ~isempty(pValid)
        pMin = min(pValid(:));
        pMax = max(pValid(:));
        if pMax > pMin
            pNormalized = (pAbs - pMin) / (pMax - pMin);
        else
            pNormalized = zeros(size(pAbs));
        end
    else
        pNormalized = zeros(size(pAbs));
    end
    pNormalized(isnan(pNormalized)) = 0;
    jetMap = colormap('jet');
    pRGB = ind2rgb(round(pNormalized * (size(jetMap, 1) - 1)) + 1, jetMap);
    
    % Blend images: pressure prevails (higher alpha for pressure)
    % Pressure transparency: 0.5 to 1.0 (more opaque for higher pressure)
    % Background shows through where pressure is low
    alphaP = pressureTransparency + (1 - pressureTransparency) * pNormalized;         % Variable transparency for pressure (0.5 to 1.0)
    
    % Blend: result = ctLayer * (1 - alphaP) + pressure * alphaP
    % When alphaP is high (high pressure), pressure dominates
    % When alphaP is low (low pressure), background shows through
    blendedRGB = ctRGB .* (1 - alphaP) + pRGB .* alphaP;
    
    figure;
    % Display blended image
    imagesc(z([1 end])*1e3, r([1 end])*1e3, blendedRGB);
    
    % Set colormap to jet for colorbar (shows pressure scale)
    colormap(gca, 'jet');
    % Add a colorbar that reflects the real pressure amplitude
    cb = colorbar;
    % Set color axis limits and ticks to match the real pressure values
    caxis([pMin pMax]);
    cb.Label.String = 'Pressure amplitude, Pa';
    % Set the colorbar colormap to jet (matches the overlay)
    colormap(cb, 'jet');
    axis equal;
    axis tight;
    xlabel(xLabelStr);
    ylabel(yLabelStr);
    title(titleStr);
    set(gca,'YDir','normal');
    
end

