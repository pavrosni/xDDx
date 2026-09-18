function plot_3d_field(xField3D, yField3D, zField3D, pField3D, pMax3D, ...
    xMax3D, yMax3D, zMax3D, levelArray, transparencyArray, toolboxLibDir, varargin)
%PLOT_3D_FIELD Plot 3D isosurfaces of pressure field.

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Add xDDx toolbox library path to reuse native visualization helpers.
if ~isempty(toolboxLibDir) && exist(toolboxLibDir, 'dir') ~= 0
    addpath(toolboxLibDir);
end

figure;
set(gcf, 'Units', 'pixels');
cmap = colormap('jet');
climits = [0 1];
scalar = levelArray;
scalarclamped = scalar;
scalarclamped(scalar < climits(1)) = climits(1);
scalarclamped(scalar > climits(2)) = climits(2);
colorsLevels = interp1(linspace(climits(1), climits(2), size(cmap, 1)), ...
    cmap, ...
    scalarclamped);

for i = 1:length(levelArray)
    p = patch(isosurface(xField3D*1e3, yField3D*1e3, zField3D*1e3, ...
        abs(pField3D)/max(abs(pField3D(:))), levelArray(i)));
    set(p, 'FaceColor', squeeze(colorsLevels(i,:)), 'EdgeColor', 'none', ...
        'FaceAlpha', transparencyArray(i));
end

% Optional transducer surface overlay (pass a TransducerSf struct).
if ~isempty(varargin)
    transducerSf = varargin{1};
    if isstruct(transducerSf) && all(isfield(transducerSf, ...
            {'xGrid', 'yGrid', 'complexVelocityAmplitude'}))
        if isfield(transducerSf, 'zGrid') && exist('check_transducer_z', 'file') == 2
            transducerSf = check_transducer_z(transducerSf);
        end

        if isvector(transducerSf.xGrid) && exist('reshape_transducer_dim', 'file') == 2
            transducerSf = reshape_transducer_dim(transducerSf);
        end

        hasStructuredGrid = ismatrix(transducerSf.xGrid) && ismatrix(transducerSf.yGrid) && ...
            ismatrix(transducerSf.zGrid) && ismatrix(transducerSf.complexVelocityAmplitude) && ...
            isequal(size(transducerSf.xGrid), size(transducerSf.yGrid), ...
            size(transducerSf.zGrid), size(transducerSf.complexVelocityAmplitude));

        if hasStructuredGrid
            cDatavAmpl = abs(transducerSf.complexVelocityAmplitude) * 1e3;
            cDataRange = max(cDatavAmpl(:)) - min(cDatavAmpl(:));
            if cDataRange > eps
                cDataNormalized = (cDatavAmpl - min(cDatavAmpl(:))) / cDataRange;
            else
                cDataNormalized = zeros(size(cDatavAmpl));
            end
            cDataNormalized = round(cDataNormalized * (size(cmap, 1) - 1)) + 1;
            rgbImageAmpl = ind2rgb(cDataNormalized, cmap);

            activeSurfaceMask = ones(size(transducerSf.complexVelocityAmplitude));
            activeSurfaceMask(abs(transducerSf.complexVelocityAmplitude) < eps) = 0;

            if ~isOctave
                hSurface = surface(transducerSf.xGrid*1e3, transducerSf.yGrid*1e3, ...
                    transducerSf.zGrid*1e3, rgbImageAmpl, ...
                    'FaceColor', 'texturemap', 'EdgeColor', 'none');
                set(hSurface, 'AlphaData', activeSurfaceMask, 'FaceAlpha', 'texturemap', ...
                    'CDataMapping', 'direct');
            else
                surface(transducerSf.xGrid*1e3, transducerSf.yGrid*1e3, ...
                    transducerSf.zGrid*1e3, rgbImageAmpl, ...
                    'FaceColor', 'texturemap', 'EdgeColor', 'none');
            end
        else
            warning(['Transducer surface is not shown because xGrid/yGrid/zGrid/complexVelocityAmplitude ', ...
                'are not consistent 2D arrays of the same size.']);
        end
    end
end

axis equal;
axis tight;
xlabel('x, mm');
ylabel('y, mm');
zlabel('z, mm');
colormap jet;
titleFor3D = ['Pressure amplitude max is ' num2str(pMax3D) ' Pa at ' ...
    '(' num2str(xMax3D*1e3) ', ' num2str(yMax3D*1e3) ', ' num2str(zMax3D*1e3) ') mm'];
cbh = colorbar;
caxis([0 1]);
set(cbh, 'XTick', flip(levelArray));
title(titleFor3D);
view(3);
if ~isOctave
    camlight;
else
    zlim([0 max(zField3D(:))*1e3]);
end

% Slice-by-slice GUI representation, as in xDDx simulation examples.
if exist('gui_plot_3d', 'file') == 2
    gui_plot_3d(xField3D, yField3D, zField3D, pField3D, levelArray, transparencyArray, false);
else
    warning('Slice-by-slice GUI is unavailable: gui_plot_3d.m was not found on MATLAB path.');
end

end
