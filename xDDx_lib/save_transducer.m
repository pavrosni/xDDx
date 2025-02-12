% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function save_transducer(isSphericalSource, expSign, frequency, xSource, ySource, zSource, radiusOfCurvature, dxSource, dySource, vSource)
TransducerSf = [];
TransducerSf.expSign = expSign;
TransducerSf.frequency = frequency;

if isSphericalSource
    TransducerSf.radiusOfCurvature = radiusOfCurvature;
end

TransducerSf.xGrid = xSource;
TransducerSf.yGrid = ySource;
TransducerSf.zGrid = zSource;
TransducerSf.dx = dxSource;
TransducerSf.dy = dySource;
TransducerSf.complexVelocityAmplitude = vSource;

currentDateTime = datetime('now', 'Format','yyyy_MM_dd_HH_mm_ss');
formattedDateTime = char(currentDateTime);

save(['transducer_' formattedDateTime '.mat'], 'TransducerSf');
end

