% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function gui_find_roc(vSource, expSign, frequency, regime, simulationDevice, isTransient, SourceParameters, FieldParameters, Medium, ServiceParameters, radiusOfCurvature, xMaxMechCoord, yMaxMechCoord, zMaxMechCoord, directionVectorMechCoord, xZeroPhase, yZeroPhase, radiusMax, dxSource, dySource)

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

xSource = SourceParameters.xGrid;
ySource = SourceParameters.yGrid;
zPrefocalShift = radiusOfCurvature - FieldParameters.zGrid(1);

TransducerSf = [];
TransducerSf.expSign = expSign;
TransducerSf.frequency = frequency;
TransducerSf.radiusOfCurvature = radiusOfCurvature;
TransducerSf.xGrid = SourceParameters.xGrid;
TransducerSf.yGrid = SourceParameters.yGrid;
TransducerSf.zGrid = SourceParameters.zGrid;
TransducerSf.dx = dxSource;
TransducerSf.dy = dySource;
TransducerSf.complexVelocityAmplitude = vSource;

lengthFigMin = 980;
widthFigMin = 450;

positionScreen = get(0,'screensize');
lengthScreen = positionScreen(3);
widthScreen = positionScreen(4);

lengthFig = max(fix(lengthScreen/2),lengthFigMin);
widthFig = max(fix(widthScreen/2), widthFigMin);
figPosition = [fix((lengthScreen-lengthFig)/2), fix((widthScreen-widthFig)/2), lengthFig, widthFig];

if isOctave
  figPosition = [100 100 lengthFigMin widthFigMin];
end

figure('Name', 'GUI: Backpropagate and Replot', 'Position', figPosition);


ax1 = subplot(1, 2, 1);
imagesc(ax1,xSource([1 end])*1e3,ySource([1 end])*1e3,abs(vSource)*1e3);
set(ax1,'YDir','normal');
colormap(ax1,'jet');
axis(ax1,'equal');

if ~isOctave
  axis(ax1,'tight');
end

title(ax1,'Vibrational velocity amplitude at the array surface, mm/s');
xlabel(ax1,'{\itx}, mm');
ylabel(ax1,'{\ity}, mm');

if isOctave
  xlim(ax1,xSource([1 end])*1e3)
  ylim(ax1,ySource([1 end])*1e3)
end



ax2 = subplot(1, 2, 2);

if (~isempty(xZeroPhase)) && (~isempty(yZeroPhase))
    [~, iZeroPhase] = min((xSource(:) - xZeroPhase).^2 + (ySource(:) - yZeroPhase).^2);

    angSource = shift_2pi(angle(vSource) - angle(vSource(iZeroPhase)));
    angSource(xSource.^2 + ySource.^2 > radiusMax^2) = 4;
else
    angSource = angle(vSource);
    angSource(xSource.^2 + ySource.^2 > radiusMax^2) = 0;
end

imagesc(ax2,xSource([1 end])*1e3,ySource([1 end])*1e3, angSource);
set(ax2,'YDir','normal');
cbh2 = colormap(ax2,'hsv');
axis(ax2,'equal');

if ~isOctave
  axis(ax2,'tight');
end

title(ax2,'Vibrational velocity phase at the array surface, rad');
xlabel(ax2,'{\itx}, mm');
ylabel(ax2,'{\ity}, mm');

if isOctave
  xlim(ax2,xSource([1 end])*1e3)
  ylim(ax2,ySource([1 end])*1e3)
end

if ~isOctave

  set (ax1,'Units', 'normalized', 'Position', [0.1 , 0.27, 0.35, 0.6]);
  set (ax2,'Units', 'normalized', 'Position', [0.6, 0.27, 0.35, 0.6]);

  uicontrol('Style', 'text', 'String', 'Radius of Curvature in mm', ...
    'Units', 'normalized', 'Position', [0.01,0.07,0.17,0.06]);

  editField = uicontrol('Style', 'edit', 'String', num2str(radiusOfCurvature*1e3), ...
    'Units', 'normalized', 'Position', [0.18,0.082525252525253,0.08,0.06]);

  uicontrol('Style', 'pushbutton', 'String', 'Backproject and Replot', ...
    'Units', 'normalized', 'Position', [0.29,0.08,0.17,0.06], ...
    'Callback', @replotCallback);

  uicontrol('Style', 'pushbutton', 'String', 'Save TransducerSf', ...
    'Units', 'normalized', 'Position', [0.47,0.08,0.17,0.06], ...
    'Callback', @saveBackpropagationCallback);

  uicontrol('Style', 'pushbutton', 'String', 'Save Alignment Parameters', ...
    'Units', 'normalized', 'Position', [0.65,0.08,0.17,0.06], ...
    'Callback', @saveRotationCallback);

  uicontrol('Style', 'pushbutton', 'String', 'Save Aligned Holo', ...
    'Units', 'normalized', 'Position', [0.83,0.08,0.11,0.06], ...
    'Callback', @saveRotatedHoloCallback);

else

  fontSizeUIControls = 8;

  uicontrol('Style', 'text', 'String', 'Radius of Curvature in mm', ...
    'Units', 'normalized', 'Position', [0.01,0.07,0.17,0.06], 'FontSize', fontSizeUIControls,'BackgroundColor','w');

  editField = uicontrol('Style', 'edit', 'String', num2str(radiusOfCurvature*1e3), ...
    'Units', 'normalized', 'Position', [0.18,0.082525252525253,0.08,0.06], 'FontSize', fontSizeUIControls,'BackgroundColor','w');

  uicontrol('Style', 'pushbutton', 'String', 'Backpropagate and Replot', ...
    'Units', 'normalized', 'Position', [0.29,0.08,0.17,0.06], 'FontSize', fontSizeUIControls, ...
    'Callback', @replotCallback);

  uicontrol('Style', 'pushbutton', 'String', 'Save Backpropagation Results', ...
    'Units', 'normalized', 'Position', [0.47,0.08,0.17,0.06], 'FontSize', fontSizeUIControls, ...
    'Callback', @saveBackpropagationCallback);

  uicontrol('Style', 'pushbutton', 'String', 'Save Alignment Parameters', ...
    'Units', 'normalized', 'Position', [0.65,0.08,0.17,0.06], 'FontSize', fontSizeUIControls, ...
    'Callback', @saveRotationCallback);

  uicontrol('Style', 'pushbutton', 'String', 'Save Rotated Holo', ...
    'Units', 'normalized', 'Position', [0.83,0.08,0.11,0.06], 'FontSize', fontSizeUIControls, ...
    'Callback', @saveRotatedHoloCallback);

end

  set (ax1,'Units', 'normalized', 'Position', [0.05 , 0.27, 0.4, 0.5]);
  set (ax2,'Units', 'normalized', 'Position', [0.5, 0.27, 0.4, 0.5]);

  colorbar(ax1);
colorbar(ax2);

    function replotCallback(~, ~)
        radiusOfCurvature = str2double(get(editField,'String'))*1e-3;

        % Check if the input is a valid number
        if ~isnan(radiusOfCurvature)
            SourceParameters.zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - SourceParameters.xGrid.^2 - SourceParameters.yGrid.^2);
            FieldParameters.zGrid = (radiusOfCurvature-zPrefocalShift)*ones(size(FieldParameters.zGrid));

            zSource = SourceParameters.zGrid;

            xSourceWithin = xSource(xSource.^2 + ySource.^2 <= radiusMax^2);
            ySourceWithin = ySource(xSource.^2 + ySource.^2 <= radiusMax^2);
            zSourceWithin = zSource(xSource.^2 + ySource.^2 <= radiusMax^2);

            %Backpropagate the rotated hologram to the surface of the source
            SourceParametersWithin.xGrid = xSourceWithin;
            SourceParametersWithin.yGrid = ySourceWithin;
            SourceParametersWithin.zGrid = zSourceWithin;


            TransducerSf.radiusOfCurvature = radiusOfCurvature;
            TransducerSf.zGrid = SourceParameters.zGrid;

            [ vSourceChangedWithin ] = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, SourceParametersWithin, FieldParameters, Medium, ServiceParameters, radiusOfCurvature);

            vSourceChanged = zeros(size(xSource));
            vSourceChanged(xSource.^2 + ySource.^2 <= radiusMax^2) = vSourceChangedWithin;

            TransducerSf.complexVelocityAmplitude = vSourceChanged;


            imagesc(ax1,xSource([1 end])*1e3,ySource([1 end])*1e3,abs(vSourceChanged)*1e3);
            set(ax1,'YDir','normal');
            colormap(ax1,'jet');
            axis(ax1,'equal');

            if ~isOctave
              axis(ax1,'tight');
            end

            title(ax1,'Vibrational velocity amplitude at the array surface, mm/s');
            xlabel(ax1,'{\itx}, mm');
            ylabel(ax1,'{\ity}, mm');

            if isOctave
              xlim(ax1,xSource([1 end])*1e3)
              ylim(ax1,ySource([1 end])*1e3)
            end

            if (~isempty(xZeroPhase)) && (~isempty(yZeroPhase))

                [~, iZeroPhase] = min((xSource(:) - xZeroPhase).^2 + (ySource(:) - yZeroPhase).^2);

                angSourceChanged = shift_2pi(angle(vSourceChanged) - angle(vSourceChanged(iZeroPhase)));
                angSourceChanged(xSource.^2 + ySource.^2 > radiusMax^2) = 4;
            else
                angSourceChanged = angle(vSourceChanged);
                angSourceChanged(xSource.^2 + ySource.^2 > radiusMax^2) = 0;
            end



            imagesc(ax2,xSource([1 end])*1e3,ySource([1 end])*1e3,angSourceChanged);
            set(ax2,'YDir','normal');
            colormap(ax2,'hsv');
            axis(ax2,'equal');

            if ~isOctave
                axis(ax2,'tight');
            end

            title(ax2,'Vibrational velocity phase at the array surface, rad');
            xlabel(ax2,'{\itx}, mm');
            ylabel(ax2,'{\ity}, mm');

            if isOctave
              xlim(ax2,xSource([1 end])*1e3)
              ylim(ax2,ySource([1 end])*1e3)
            end

            set (ax1,'Units', 'normalized', 'Position', [0.05 , 0.27, 0.4, 0.5]);
            set (ax2,'Units', 'normalized', 'Position', [0.5, 0.27, 0.4, 0.5]);

            cbh1 = colorbar(ax1);
            cbh2 = colorbar(ax2);

        else
            errordlg('Invalid input. Please enter a valid number.', 'Error');
        end
    end


    function saveBackpropagationCallback(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0
            save(fullfile(path, filename), 'TransducerSf');
            disp(['Variable BackpropagationResults saved as ', fullfile(path, filename)]);
        end
    end

    function saveRotationCallback(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0
            zMaxAcoustCoord = radiusOfCurvature;
            save(fullfile(path, filename), 'zMaxAcoustCoord', 'xMaxMechCoord','yMaxMechCoord','zMaxMechCoord', 'directionVectorMechCoord');
            disp(['Hologram position parameters saved as ', fullfile(path, filename)]);
        end
    end

    function saveRotatedHoloCallback(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0

            HologramSf = [];
            HologramSf.expSign = expSign;
            HologramSf.frequency = frequency;
            HologramSf.xGrid = FieldParameters.xGrid;
            HologramSf.yGrid = FieldParameters.yGrid;
            HologramSf.zPosition = radiusOfCurvature-zPrefocalShift;
            HologramSf.dx = FieldParameters.dx;
            HologramSf.dy = FieldParameters.dy;
            HologramSf.complexPressureAmplitude = FieldParameters.input;

            save(fullfile(path, filename), 'HologramSf', 'Medium', 'radiusOfCurvature');
            disp(['Rotated hologram saved as ', fullfile(path, filename)]);
        end
    end

end
