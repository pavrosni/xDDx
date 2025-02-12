% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function gui_plot_3d(xField3D,yField3D, zField3D, pField3D, levelArray, transparencyArray, show3D)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if nargin < 7
    show3D = true;
end

%Find the pressure maximum point (xMax3D, yMax3D, zMax3D)
[pMax3D,iMax3D] = max(abs(pField3D(:)));
xMax3D = xField3D(iMax3D);
yMax3D = yField3D(iMax3D);
zMax3D = zField3D(iMax3D);

xFieldVector = squeeze(xField3D(1,:,1));

%Plot the 3D field isosurfaces given by 'levelArray' with transparencies given by 'transparencyArray'
% positionCurrent = get(groot,'DefaultFigurePosition');
if show3D
    figure;
    hold on;
    cmap = colormap('jet');
    climits = [0 1];
    scalar=levelArray;
    scalarclamped = scalar;
    scalarclamped(scalar < climits(1)) = climits(1);
    scalarclamped(scalar > climits(2)) = climits(2);
    colorsLevels = interp1(linspace(climits(1), climits(2), size(cmap, 1)), ...
        cmap, ...
        scalarclamped);

    for i = 1:length(levelArray)
        p=patch(isosurface(xField3D*1e3, yField3D*1e3, zField3D*1e3, abs(pField3D)/max(abs(pField3D(:))), levelArray(i)));
        set(p,'FaceColor',squeeze(colorsLevels(i,:)),'EdgeColor','none','FaceAlpha',transparencyArray(i));
    end
    axis equal;
    axis tight;
    xlabel('x, mm');
    ylabel('y, mm');
    zlabel('z, mm');
    title(['Pressure amplitude max is ' num2str(pMax3D) ' Pa at ' '(' num2str(xMax3D*1e3) ', ' num2str(yMax3D*1e3) ', ' num2str(zMax3D*1e3) ') mm' ]);
    colormap jet;

    titleFor3D = ['Pressure amplitude max is ' num2str(pMax3D) ' Pa at ' '(' num2str(xMax3D*1e3) ', ' num2str(yMax3D*1e3) ', ' num2str(zMax3D*1e3) ') mm' ];

    if ~isOctave
        cbh = colorbar;
        caxis([0 1]);
        set(cbh,'XTick',flip(levelArray));
        view(3);
        camlight;
    else
        titleFor3D = {titleFor3D, ['pressure amplitude iso-levels: ' num2str(levelArray)]};
        view(3);
    end

    title(titleFor3D);


end

%Plot the 2D slices passing through the zy and zx axial planes
[~, ijkMax] = max(abs(pField3D(:)));
[~, jMax, ~] = ind2sub(size(pField3D),ijkMax);
xField2D = squeeze(xField3D(:,jMax,:));
yField2D = squeeze(yField3D(:,jMax,:));
zField2D = squeeze(zField3D(:,jMax,:));
pField2D = squeeze(pField3D(:,jMax,:));
xCurrent = xField3D(ijkMax);
jCurrent = jMax;

positionScreen = get(0,'screensize');
lengthScreen = positionScreen(3);
widthScreen = positionScreen(4);

lengthFig = fix(lengthScreen/2);
widththFig = fix(widthScreen/2);
figPosition = [fix((lengthScreen-lengthFig)/2), fix((widthScreen-widththFig)/2), lengthFig, widththFig];

if isOctave
    figPosition = [100 100 720 450];
end

figure('Name', 'GUI: 3D Field', 'Position', figPosition);



labelX = uicontrol('Style', 'text', 'String', 'x-position in mm', ...
    'Units', 'normalized', 'Position', [0.01,0.02,0.17,0.06]);

editField = uicontrol('Style', 'edit', 'String', num2str(xCurrent*1e3), ...
    'Units', 'normalized', 'Position', [0.18,0.032525252525253,0.08,0.06]);

uicontrol('Style', 'pushbutton', 'String', 'Replot', ...
    'Units', 'normalized', 'Position', [0.29,0.03,0.17,0.06], ...
    'Callback', @replotCallback);

labelChangeX = uicontrol('Style', 'text', 'String', 'Change x', ...
    'Units', 'normalized', 'Position', [0.49,0.068,0.17,0.031]);


slider = uicontrol('Style', 'slider', ...
    'Units', 'normalized',...
    'Min', 1, 'Max', length(xFieldVector), 'Value', jMax, ...
    'SliderStep',[1/(length(xFieldVector)-1),1/(length(xFieldVector)-1)], ...
    'Position', [0.47,0.03,0.2,0.03], ...
    'Callback', @sliderCallback);

uicontrol('Style', 'pushbutton', 'String', 'Save Layer', ...
    'Units', 'normalized', 'Position', [0.69,0.03,0.146038338658149,0.06], ...
    'Callback', @saveLayer);

uicontrol('Style', 'pushbutton', 'String', 'Save 3D', ...
    'Units', 'normalized', 'Position', [0.85,0.03,0.11,0.06], ...
    'Callback', @save3D);

if ~isOctave
    addlistener(slider, 'ContinuousValueChange', @sliderCallback);
end

if isOctave
    set(labelX,'BackgroundColor','w');
    set(editField,'BackgroundColor','w');
    set(labelChangeX,'BackgroundColor','w');
end

ax = axes;
replotZy;


    function replotCallback(~, ~)



        xCurrent = str2double(get(editField,'String'))*1e-3;
        [~, jCurrent] = min(abs(xFieldVector - xCurrent));

        if isnan(xCurrent)
            errordlg('Invalid input. Please enter a valid number.', 'Error');
            return;
        end

        replotZy
    end


    function saveLayer(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0
            save(fullfile(path, filename), 'xField2D', 'yField2D','zField2D','pField2D');
            disp(['2D Field saved as ', fullfile(path, filename)]);
        end
    end


    function save3D(~, ~)
        [filename, path] = uiputfile('*.mat', 'Save as');
        if filename ~= 0
            save(fullfile(path, filename), 'xField3D', 'yField3D','zField3D','pField3D');
            disp(['3D Field saved as ', fullfile(path, filename)]);
        end
    end

    function replotZy

        xField2D = squeeze(xField3D(:,jCurrent,:));
        yField2D = squeeze(yField3D(:,jCurrent,:));
        zField2D = squeeze(zField3D(:,jCurrent,:));
        pField2D = squeeze(pField3D(:,jCurrent,:));

        xCurrent = xFieldVector(jCurrent);

        set(editField,'String',num2str(xCurrent*1e3));
        set(slider, 'Value', jCurrent);

        if ~isOctave
            contourf(ax, zField2D*1e3, yField2D*1e3, abs(pField2D), 100, 'LineStyle', 'none');
        else
            imagesc(ax, zField2D([1 end])*1e3, yField2D([1 end])*1e3, abs(pField2D));
            set(ax,'YDir','normal');
        end


        title(ax, ['Pressure amplitude zy-distribution in Pa (x = ' num2str(xField2D(1)*1e3) ' mm)' ]);
        xlabel(ax, 'z, mm');
        ylabel(ax, 'y, mm');
        colormap(ax, 'jet');
        axis(ax, 'equal');

        positionAxes = get(ax,'Position');
        axesRatioLonH = positionAxes(3)/positionAxes(4);

        if ~isOctave
            heightsAxes = 0.7550;
            posxAxes = 0.1229;
            posyAxes = 0.19;
        else
            heightsAxes = 0.6550;
            posxAxes = 0.2300;
            posyAxes = 0.2100;
        end
        set(ax,'Position',[posxAxes    posyAxes    axesRatioLonH*heightsAxes    heightsAxes]);


        if isOctave
            xlim(ax,zField2D([1 end])*1e3)
            ylim(ax,yField2D([1 end])*1e3)
        end

        colorbar(ax);

        %        if ~isOctave
        %         set(ax,'Position',[0.1229    0.19    0.7325    0.7550]);
        %        else
        %         set(ax,'Position',[0.1300   0.1100   0.6200   0.8150]);
        %      end
    end

    function sliderCallback(src, ~)

        jCurrent = round(get(src, 'Value'));

        replotZy

    end
end
