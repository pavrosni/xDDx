% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xField2D, yField2D, zField2D, pField2D] = plot_2d_field(xField3D,yField3D, zField3D, pField3D)
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

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

figure;
if ~isOctave
    contourf(outX*1e3, outY*1e3, abs(pField2D), 100, 'LineStyle', 'none');
else
    imagesc(outX([1 end])*1e3, outY([1 end])*1e3, abs(pField2D));
    set(gca,'YDir','normal');
end

title(['Pressure amplitude ' outputAxes '-distribution in Pa (' constAxis ' = ' num2str(constAxisValue*1e3) ' mm)' ]);
xlabel([outputAxes(1) ', mm']);
ylabel([outputAxes(2) ', mm']);
colormap('jet');
colorbar;
axis('equal');

end

