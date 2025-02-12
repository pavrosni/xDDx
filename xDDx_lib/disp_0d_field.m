% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [xField0D, yField0D, zField0D, pField0D] = disp_0d_field(xField3D,yField3D, zField3D, pField3D)

xField0D = squeeze(xField3D);
yField0D = squeeze(yField3D);
zField0D = squeeze(zField3D);
pField0D = squeeze(pField3D);


disp(['Pressure amplitude at x = ' num2str(xField0D*1e3) ' mm, y = ' num2str(yField0D*1e3) ' mm, z = ' num2str(zField0D*1e3) ' mm is ' num2str(abs(pField0D)) ' Pa']);

end

