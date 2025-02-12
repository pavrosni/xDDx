% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function gui_transient_frame_by_frame(xSource, ySource, vSourceFrames, frameStep, indexFirstFrame, indexLastFrame, isTimeDomain, radiusOfCurvature, dxSource, dySource,  expSign, cascadeShift)

isSphericalSource = true;
if isempty(radiusOfCurvature)
    isSphericalSource = false;
end

if nargin < 12
    cascadeShift = 0;
end

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
idxFrameCurrent = indexFirstFrame;

if isTimeDomain
    transientType = 'Time Domain';
    variableType = 't';
    titleMainPart = 'Vibrational velocity [mm/s] at {\itt} =';
    titlePostfix = ' us';
    buttonTextSave2D = 'Save t-sample';
    buttonTextSave3D = 'Save TransducerTr';

    vFrameCurrent = squeeze(vSourceFrames(:,:,idxFrameCurrent));
    frameCurrent = (idxFrameCurrent-indexFirstFrame)*frameStep;
    titleFrameInfo = num2str(frameCurrent*1e6, '%1.4f');
    editFieldString = ['Time in ' titlePostfix];
    minValue = min(vSourceFrames(:));
    maxValue = max(vSourceFrames(:));
else
    transientType = 'Frequency Domain';
    variableType = 'f';
    buttonTextSave2D = 'Save  TransducerSf';
    buttonTextSave3D = 'Save all TransducerSf';

    titleMainPart = 'Vibrational velocity amplitude [mm/s] at {\itf} =';
    titlePostfix = ' MHz';
    vFrameCurrent = squeeze(abs(vSourceFrames(:,:,idxFrameCurrent)));
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


    uicontrol('Style', 'pushbutton', 'String', buttonTextSave2D, ...
        'Units', 'normalized', 'Position', [0.69,0.03,0.146038338658149,0.06], ...
        'Callback', @saveLayer);

    uicontrol('Style', 'pushbutton', 'String', buttonTextSave3D, ...
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

    uicontrol('Style', 'pushbutton', 'String', buttonTextSave2D, ...
        'Units', 'normalized', 'Position', [0.69,0.03,0.146038338658149,0.06], ...
        'Callback', @saveLayer, 'FontSize', fontSizeUIControls);

    uicontrol('Style', 'pushbutton', 'String', buttonTextSave3D, ...
        'Units', 'normalized', 'Position', [0.85,0.03,0.11,0.06], ...
        'Callback', @save3D, 'FontSize', fontSizeUIControls);


end

ax = axes;

replotFrame;

    function replotCallback(~, ~)

        if isTimeDomain
            idxFrameCurrent = round(str2double(get(editField,'String'))*1e-6/frameStep) + indexFirstFrame;
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
                SurfaceVelocitySample = [];
                SurfaceVelocitySample.time = (idxFrameCurrent-indexFirstFrame)*frameStep;

                if isSphericalSource
                    SurfaceVelocitySample.radiusOfCurvature = radiusOfCurvature;
                end

                SurfaceVelocitySample.xGrid = xSource;
                SurfaceVelocitySample.yGrid = ySource;

                if isSphericalSource
                    SurfaceVelocitySample.zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xSource.^2 - ySource.^2);
                else
                    SurfaceVelocitySample.zGrid = zeros(size(xSource));
                end

                SurfaceVelocitySample.dx = dxSource;
                SurfaceVelocitySample.dy = dySource;                
                SurfaceVelocitySample.velocity = vSourceFrames(:,:,idxFrameCurrent);
                save(fullfile(path, filename), 'SurfaceVelocitySample');
            else
                TransducerSf = [];
                TransducerSf.expSign = expSign;
                TransducerSf.frequency = (idxFrameCurrent - indexFirstFrame + 1)*frameStep;

                if isSphericalSource
                    TransducerSf.radiusOfCurvature = radiusOfCurvature;
                end

                TransducerSf.xGrid = xSource;
                TransducerSf.yGrid = ySource;

                if isSphericalSource
                    TransducerSf.zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xSource.^2 - ySource.^2);
                else
                    TransducerSf.zGrid = zeros(size(xSource));
                end

                TransducerSf.dx = dxSource;
                TransducerSf.dy = dySource;
                TransducerSf.complexVelocityAmplitude = squeeze(vSourceFrames(:,:,idxFrameCurrent));
                save(fullfile(path, filename), 'TransducerSf');
            end

            disp(['Data saved as ', fullfile(path, filename)]);
        end
    end


    function save3D(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0

            if isTimeDomain
                TransducerTr = [];
                TransducerTr.time = (0:(indexLastFrame - indexFirstFrame))*frameStep;

                if isSphericalSource
                    TransducerTr.radiusOfCurvature = radiusOfCurvature;
                end

                TransducerTr.xGrid = xSource;
                TransducerTr.yGrid = ySource;

                if isSphericalSource
                    TransducerTr.zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xSource.^2 - ySource.^2);
                else
                    TransducerTr.zGrid = zeros(size(xSource));
                end

                TransducerTr.dx = dxSource;
                TransducerTr.dy = dySource;
                TransducerTr.velocity = vSourceFrames(:,:,indexFirstFrame:indexLastFrame);
                save(fullfile(path, filename), 'TransducerTr', '-v7.3');
            else
                TransducersSf = [];

                for iFreq = indexFirstFrame:indexLastFrame
                    iTx = iFreq - indexFirstFrame + 1;
                    TransducersSf{iTx}.frequency = iFreq*frameStep;

                    if isSphericalSource
                        TransducersSf{iTx}.radiusOfCurvature = radiusOfCurvature;
                    end

                    TransducersSf{iTx}.xGrid = xSource;
                    TransducersSf{iTx}.yGrid = ySource;

                    if isSphericalSource
                        TransducersSf{iTx}.zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - xSource.^2 - ySource.^2);
                    else
                        TransducersSf{iTx}.zGrid = zeros(size(xSource));
                    end

                    TransducersSf{iTx}.dx = dxSource;
                    TransducersSf{iTx}.dy = dySource;
                    TransducersSf{iTx}.complexVelocityAmplitudes = vSourceFrames(:,:,iFreq);
                end

                save(fullfile(path, filename), 'TransducersSf', '-v7.3');
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
            frameCurrent = (idxFrameCurrent-indexFirstFrame)*frameStep;
            titleFrameInfo = num2str(frameCurrent*1e6, '%1.4f');
            vFrameCurrent = squeeze(vSourceFrames(:,:,idxFrameCurrent));
        else
            frameCurrent = idxFrameCurrent*frameStep;
            titleFrameInfo = num2str(frameCurrent*1e-6, '%1.4f');
            vFrameCurrent = squeeze(abs(vSourceFrames(:,:,idxFrameCurrent)));
        end

        set(editField,'String',num2str(titleFrameInfo));

        imagesc(ax, xSource([1 end])*1e3, ySource([1 end])*1e3, vFrameCurrent*1e3);
        set(ax,'YDir','normal');

        title(ax, [titleMainPart titleFrameInfo titlePostfix]);
        xlabel(ax, 'z, mm');
        ylabel(ax, 'y, mm');

        if isTimeDomain
            caxis([minValue maxValue]*1e3);
        end

        colormap(ax, 'jet');
        axis(ax, 'equal');
        axis(ax, 'tight');

        if isOctave
            xlim(ax,xSource([1 end])*1e3)
            ylim(ax,ySource([1 end])*1e3)
        end

        set(ax, 'Position', [0.1425    0.2100    0.7120    0.7150]);

        colorbar(ax);

    end

    function sliderCallback(src, ~)

        idxFrameCurrent = round(get(src, 'Value'));

        replotFrame

    end


end
