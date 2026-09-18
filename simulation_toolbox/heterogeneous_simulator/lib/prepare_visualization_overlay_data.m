function [MaterialMatrixForPlot, xGridMediumVecForPlot, yGridMediumVecForPlot, zGridMediumVecForPlot] = ...
    prepare_visualization_overlay_data(MaterialMatrix, xGridMediumVec, yGridMediumVec, zGridMediumVec, strongMemorySavingMode)
%PREPARE_VISUALIZATION_OVERLAY_DATA Keep only the medium data needed for plots.

if strongMemorySavingMode
    MaterialMatrixForPlot = [];
    xGridMediumVecForPlot = [];
    yGridMediumVecForPlot = [];
    zGridMediumVecForPlot = [];
else
    MaterialMatrixForPlot = struct('c0', MaterialMatrix.c0);
    xGridMediumVecForPlot = xGridMediumVec;
    yGridMediumVecForPlot = yGridMediumVec;
    zGridMediumVecForPlot = zGridMediumVec;
end

end
