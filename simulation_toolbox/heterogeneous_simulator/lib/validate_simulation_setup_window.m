function proceed = validate_simulation_setup_window(MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, aperture, radiusOfCurvature, radialReserveX, radialReserveY, xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, xGridSimulationVec, yGridSimulationVec, zGridSimulationVec)
%VALIDATE_SIMULATION_SETUP_WINDOW Show yz/xz sketches for user validation.
%
% radialReserveX, radialReserveY: reserve factors for the boundary transverse sizes.
% Returns:
%   proceed (logical): true if user clicks Proceed, false otherwise.

proceed = false;

fig = figure( ...
    'Name', 'Validate simulation setup', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'figure', ...
    'Units', 'normalized', ...
    'Position', [0.15 0.15 0.7 0.7], ...
    'Color', 'w', ...
    'CloseRequestFcn', @onClose);

% Reserve room below each axes for colorbar + legend and keep controls clear.
ax1 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.08 0.29 0.38 0.65]);
ax2 = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.56 0.29 0.38 0.65]);

plot_simulation_setup_sketch(ax1, 'yz', MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, ...
    aperture, radiusOfCurvature, xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, ...
    xGridSimulationVec, yGridSimulationVec, zGridSimulationVec, true, radialReserveX, radialReserveY, true, true);
plot_simulation_setup_sketch(ax2, 'xz', MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, ...
    aperture, radiusOfCurvature, xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, ...
    xGridSimulationVec, yGridSimulationVec, zGridSimulationVec, true, radialReserveX, radialReserveY, false);

uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.05 0.00 0.55 0.04], ...
    'String', 'Validate target, boundary condition, and simulation region. Closing window counts as Stop.', ...
    'HorizontalAlignment', 'left', 'BackgroundColor', 'w');

hKeepOpen = uicontrol('Parent', fig, 'Style', 'checkbox', 'Units', 'normalized', ...
    'Position', [0.05 0.045 0.40 0.03], ...
    'String', 'Keep window open after Proceed/Stop', ...
    'Value', 0, ...
    'BackgroundColor', 'w', ...
    'Callback', @onToggleKeepOpen);

setappdata(fig, 'kw_keep_open', false);

hStop = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.64 0.015 0.15 0.045], ...
    'String', 'Stop', 'FontWeight', 'bold', ...
    'Callback', @onStop);

hProceed = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.81 0.015 0.15 0.045], ...
    'String', 'Proceed', 'FontWeight', 'bold', ...
    'Callback', @onProceed);

drawnow;
figure(fig);  % Bring window to front (above other windows)

try
    uiwait(fig);
catch
    proceed = false;
end

if isgraphics(fig, 'figure')
    keepOpen = getappdata(fig, 'kw_keep_open');
    if ~keepOpen
        delete(fig);
    end
end

    function onProceed(~, ~)
        if ~isgraphics(fig, 'figure')
            return;
        end
        proceed = true;
        finalizeChoiceAndMaybeClose();
    end

    function onStop(~, ~)
        if ~isgraphics(fig, 'figure')
            return;
        end
        proceed = false;
        finalizeChoiceAndMaybeClose();
    end

    function onClose(~, ~)
        if ~isgraphics(fig, 'figure')
            return;
        end
        proceed = false;
        uiresume(fig);
        delete(fig);
    end

    function onToggleKeepOpen(~, ~)
        if ~isgraphics(fig, 'figure')
            return;
        end
        setappdata(fig, 'kw_keep_open', logical(get(hKeepOpen, 'Value')));
    end

    function finalizeChoiceAndMaybeClose()
        if ~isgraphics(fig, 'figure')
            return;
        end
        keepOpen = getappdata(fig, 'kw_keep_open');
        if keepOpen
            set(hStop, 'Enable', 'off');
            set(hProceed, 'Enable', 'off');
        end
        uiresume(fig);
        if ~keepOpen && isgraphics(fig, 'figure')
            delete(fig);
        end
    end

end
