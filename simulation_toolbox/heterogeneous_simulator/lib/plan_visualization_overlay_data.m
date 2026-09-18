function [keepOverlayData, clearMediumDataEarly] = plan_visualization_overlay_data(strongMemorySavingMode, useGUI)
%PLAN_VISUALIZATION_OVERLAY_DATA Decide when medium data can be discarded.

keepOverlayData = ~strongMemorySavingMode;
clearMediumDataEarly = strongMemorySavingMode && ~useGUI;

end
