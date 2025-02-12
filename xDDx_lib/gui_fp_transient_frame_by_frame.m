% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [outX, outY, pField2D, outputAxes, constAxis, constAxisValue] = gui_fp_transient_frame_by_frame(xField3D,yField3D, zField3D, pField3D, frameStep, indexFirstFrame, indexLastFrame, isTimeDomain, radiusOfCurvature, cascadeShift)

xField2D = squeeze(xField3D);
yField2D = squeeze(yField3D);
zField2D = squeeze(zField3D);
pField2D = squeeze(pField3D);

outputAxes = [];
constAxis = [];
constAxisValue = [];
outX = [];
outY = [];
if abs(range(xField2D(:))) < eps
    outputAxes = 'zy';
    constAxis = 'x';
    constAxisValue = xField2D(1);
    outX = zField2D;
    outY = yField2D;
elseif abs(range(yField2D(:))) < eps
    outputAxes = 'zx';
    constAxis = 'y';
    constAxisValue = yField2D(1);
    outX = zField2D;
    outY = xField2D;
elseif abs(range(zField2D(:))) < eps
    outputAxes = 'xy';
    constAxis = 'z';
    constAxisValue = zField2D(1);
    outX = xField2D;
    outY = yField2D;
end


isSphericalSource = true;
if isempty(radiusOfCurvature)
    isSphericalSource = false;
end

if nargin < 10
    cascadeShift = 0;
end

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
idxFrameCurrent = indexFirstFrame;

if isTimeDomain
    transientType = 'Time Domain';
    variableType = 't';
    % titleMainPart = 'Pressure [mm/s] at {\itt} =';
    titleMainPart = ['Pressure ' outputAxes '-distribution in Pa (' constAxis ' = ' num2str(constAxisValue*1e3) ' mm) at {\itt} ='];
    titlePostfix = ' us';
    pFrameCurrent = squeeze(pField2D(:,:,idxFrameCurrent));
    frameCurrent = (idxFrameCurrent-1)*frameStep;
    titleFrameInfo = num2str(frameCurrent*1e6, '%1.4f');
    editFieldString = ['Time in ' titlePostfix];
    minValue = min(pField2D(:));
    maxValue = max(pField2D(:));
else
    transientType = 'Frequency Domain';
    variableType = 'f';
    titleMainPart = 'Pressure amplitude [mm/s] at {\itf} =';
    titlePostfix = ' MHz';
    pFrameCurrent = squeeze(abs(pField2D(:,:,idxFrameCurrent)));
    frameCurrent = idxFrameCurrent*frameStep;
    titleFrameInfo = num2str(frameCurrent*1e-6, '%1.4f');
    editFieldString = ['Frequency in ' titlePostfix];
end


positionScreen = get(0,'screensize');
lengthScreen = positionScreen(3);
widthScreen = positionScreen(4);

lengthFig = fix(lengthScreen/2);
widththFig = fix(widthScreen/2);
figPosition = [fix((lengthScreen-lengthFig)/2)+cascadeShift, fix((widthScreen-widththFig)/2)+cascadeShift, lengthFig, widththFig];

if isOctave
  figPosition = [100 100 720 450];   
end

figure('Name', ['GUI: Transient ' transientType], 'Position', figPosition);


if ~isOctave

    uicontrol('Style', 'text', 'String', editFieldString, ...
        'Units', 'normalized', 'Position', [0.01,0.02,0.17,0.06]);

    editField = uicontrol('Style', 'edit', 'String', num2str(titleFrameInfo), ...
        'Units', 'normalized', 'Position', [0.18,0.032525252525253,0.08,0.06]);

    uicontrol('Style', 'pushbutton', 'String', 'Replot', ...
        'Units', 'normalized', 'Position', [0.29,0.03,0.17,0.06], ...
        'Callback', @replotCallback);

    uicontrol('Style', 'text', 'String', ['Change ' variableType], ...
        'Units', 'normalized', 'Position', [0.49,0.068,0.17,0.031]);


    slider = uicontrol('Style', 'slider', ...
        'Units', 'normalized',...
        'Min', indexFirstFrame, 'Max', indexLastFrame, 'Value', indexFirstFrame, ...
        'SliderStep',[1/(indexLastFrame-indexFirstFrame),1/(indexLastFrame-indexFirstFrame)], ...
        'Position', [0.47,0.03,0.2,0.03], ...
        'Callback', @sliderCallback);


    uicontrol('Style', 'pushbutton', 'String', ['Save current ' variableType '-sample'], ...
        'Units', 'normalized', 'Position', [0.69,0.03,0.146038338658149,0.06], ...
        'Callback', @saveLayer);

    uicontrol('Style', 'pushbutton', 'String', ['Save all ' variableType '-samples'], ...
        'Units', 'normalized', 'Position', [0.85,0.03,0.11,0.06], ...
        'Callback', @save3D);

    addlistener(slider, 'ContinuousValueChange', @sliderCallback);


else

    fontSizeUIControls = 8;


    uicontrol('Style', 'text', 'String', editFieldString, ...
        'Units', 'normalized', 'Position', [0.01,0.02,0.17,0.06], 'FontSize', fontSizeUIControls,'BackgroundColor','w');

    editField = uicontrol('Style', 'edit', 'String', num2str(titleFrameInfo), ...
        'Units', 'normalized', 'Position', [0.18,0.032525252525253,0.08,0.06], 'FontSize', fontSizeUIControls,'BackgroundColor','w');

    uicontrol('Style', 'pushbutton', 'String', 'Replot', ...
        'Units', 'normalized', 'Position', [0.29,0.03,0.17,0.06], ...
        'Callback', @replotCallback, 'FontSize', fontSizeUIControls);

    uicontrol('Style', 'text', 'String', 'Change frame', ...
        'Units', 'normalized', 'Position', [0.49,0.068,0.17,0.021], 'FontSize', fontSizeUIControls,'BackgroundColor','w');

    slider = uicontrol('Style', 'slider', ...
        'Units', 'normalized',...
        'Min', indexFirstFrame, 'Max', indexLastFrame, 'Value', indexFirstFrame, ...
        'SliderStep',[1/(indexLastFrame-indexFirstFrame),1/(indexLastFrame-indexFirstFrame)], ...
        'Position', [0.47,0.03,0.2,0.03], ...
        'Callback', @sliderCallback, 'FontSize', fontSizeUIControls);

    uicontrol('Style', 'pushbutton', 'String', ['Save current ' variableType '-point'], ...
        'Units', 'normalized', 'Position', [0.69,0.03,0.146038338658149,0.06], ...
        'Callback', @saveLayer, 'FontSize', fontSizeUIControls);

    uicontrol('Style', 'pushbutton', 'String', ['Save all ' variableType '-points'], ...
        'Units', 'normalized', 'Position', [0.85,0.03,0.11,0.06], ...
        'Callback', @save3D, 'FontSize', fontSizeUIControls);


end

ax = axes;

replotFrame;

    function replotCallback(~, ~)

        if isTimeDomain
            idxFrameCurrent = round(str2double(get(editField,'String'))*1e-6/frameStep) + 1;
        else
            idxFrameCurrent = round(str2double(get(editField,'String'))*1e6/frameStep);
        end

        if isnan(idxFrameCurrent)
            errordlg('Invalid input. Please enter a valid number.', 'Error');
            return;
        end

        replotFrame
    end

    function saveLayer(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0
            if isTimeDomain
                % SurfaceVelocitySample = [];
                timeFrame = (idxFrameCurrent-1)*frameStep;
                pFieldFrame2D = pField2D(:,:,idxFrameCurrent);
                save(fullfile(path, filename), 'xField2D', 'yField2D', 'zField2D', 'timeFrame', 'pFieldFrame2D');
            else

            end

            disp(['Data saved as ', fullfile(path, filename)]);
        end
    end


    function save3D(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0

            if isTimeDomain
                time = ((indexFirstFrame:indexLastFrame)-1)*frameStep;
                pFieldFrames3D = pField3D(:,:,:,indexFirstFrame:indexLastFrame);
                save(fullfile(path, filename), 'xField3D', 'yField3D', 'zField3D', 'time', 'pFieldFrames3D', '-v7.3');
            else

            end

            disp(['Data saved as ', fullfile(path, filename)]);
        end
    end

    function replotFrame
        if idxFrameCurrent > indexLastFrame
            idxFrameCurrent = indexLastFrame;
        end

        if idxFrameCurrent < indexFirstFrame
            idxFrameCurrent = indexFirstFrame;
        end

        set(slider, 'Value', idxFrameCurrent);


        if isTimeDomain
            frameCurrent = (idxFrameCurrent-1)*frameStep;
            titleFrameInfo = num2str(frameCurrent*1e6, '%1.4f');
            pFrameCurrent = squeeze(pField2D(:,:,idxFrameCurrent));
        else
            frameCurrent = idxFrameCurrent*frameStep;
            titleFrameInfo = num2str(frameCurrent*1e-6, '%1.4f');
            pFrameCurrent = squeeze(abs(pField2D(:,:,idxFrameCurrent)));
        end

        set(editField,'String',num2str(titleFrameInfo));

        imagesc(ax, outX([1 end])*1e3, outY([1 end])*1e3, pFrameCurrent);
        set(ax,'YDir','normal');

        title(ax, [titleMainPart titleFrameInfo titlePostfix]);
        xlabel([outputAxes(1) ', mm']);
        ylabel([outputAxes(2) ', mm']);

        if isTimeDomain
            caxis([minValue maxValue]);
        end

        colormap(ax, 'jet');
        axis(ax, 'equal');
        axis(ax, 'tight');

        if isOctave
            xlim(ax,outX([1 end])*1e3)
            ylim(ax,outY([1 end])*1e3)
        end

        set(ax, 'Position', [0.1425    0.2100    0.7120    0.7150]);

        colorbar(ax);

    end

    function sliderCallback(src, ~)

        idxFrameCurrent = round(get(src, 'Value'));

        replotFrame

    end


end
