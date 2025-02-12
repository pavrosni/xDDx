% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xField1D, yField1D, zField1D, pField1D] = plot_1d_field(xField3D,yField3D, zField3D, pField3D)

xField1D = squeeze(xField3D);
yField1D = squeeze(yField3D);
zField1D = squeeze(zField3D);
pField1D = squeeze(pField3D);

outputAxis = [];
outX = [];
if abs(range(xField1D(:))) > eps
    outputAxis = 'x';
    outX = xField1D;
elseif abs(range(yField1D(:))) > eps
    outputAxis = 'y';
    outX = yField1D;
elseif abs(range(zField1D(:))) > eps
    outputAxis = 'z';
    outX = zField1D;
end

figure;
plot(outX*1e3, abs(pField1D), 'k');

title(['Pressure amplitude ' outputAxis '-distribution ']);
xlabel([outputAxis ', mm']);
ylabel('Pressure amplitude, Pa');
grid on;
grid minor;

end

