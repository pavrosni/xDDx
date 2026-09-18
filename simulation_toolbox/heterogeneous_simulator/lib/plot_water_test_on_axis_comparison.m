function plot_water_test_on_axis_comparison(complexPressure, ...
    xGridSimulationVec, yGridSimulationVec, zGridSimulationVec, ...
    xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, ...
    TransducerSf, soundSpeedContactMedium, densityContactMedium, ...
    frequency, aperture, isSphericalSource, radiusOfCurvature, ~, ~)
%PLOT_WATER_TEST_ON_AXIS_COMPARISON  k-Wave vs. O'Neil on-axis pressure along z (x≈0, y≈0).
%   On-axis analytical formulas match flat_piston_simulation_sf.m and
%   spherical_piston_simulation_sf.m (validation toolbox).

ixS = find((single(xFieldBegin) <= single(xGridSimulationVec)) & (single(xGridSimulationVec) <= single(xFieldEnd)));
iyS = find((single(yFieldBegin) <= single(yGridSimulationVec)) & (single(yGridSimulationVec) <= single(yFieldEnd)));
izS = find((single(zFieldBegin) <= single(zGridSimulationVec)) & (single(zGridSimulationVec) <= single(zFieldEnd)));
if isempty(ixS) || isempty(iyS) || isempty(izS)
    error('plot_water_test_on_axis_comparison: no grid points in the field box along x, y, or z.');
end

lix = find(single(xGridSimulationVec(ixS)) == single(0), 1, 'first');
liy = find(single(yGridSimulationVec(iyS)) == single(0), 1, 'first');
if isempty(lix) || isempty(liy)
    warning('plot_water_test_on_axis_comparison:noAxisPoint', ...
        ['Skipping water-test on-axis comparison because the simulation field box does not contain ' ...
         'a grid point exactly at x = 0 and y = 0 in single precision.']);
    return;
end

pKW = squeeze(complexPressure(lix, liy, :));
pKW = pKW(:);
zPlot = zGridSimulationVec(izS);
zPlot = zPlot(:);

v0 = max(abs(TransducerSf.complexVelocityAmplitude(:)));
p0 = densityContactMedium * soundSpeedContactMedium * v0;
k = 2 * pi * frequency / soundSpeedContactMedium;
a = aperture / 2;
expSign = 1;
if isfield(TransducerSf, 'expSign')
    expSign = TransducerSf.expSign;
end

zAx = zPlot;

if isSphericalSource
    F = radiusOfCurvature;
    Rmax = F * (1 + (1 - zAx / F).^2 - 2 * (1 - zAx / F) * sqrt(1 - (a / F)^2)).^0.5;
    pAnalytical = p0 ./ (1 - zAx / F) .* (exp(-1i * k * zAx) - exp(-1i * k * Rmax));
    pAnalytical(abs(zAx - F) < eps) = 1i * p0 * k * F * (1 - sqrt(1 - (a / F)^2)) * exp(-1i * k * F);
else
    pAnalytical = 2 * 1i * p0 * exp(-1i * k / 2 * (sqrt(a^2 + zAx.^2) + zAx)) .* sin(k / 2 * (sqrt(a^2 + zAx.^2) - zAx));
end
pAnalytical = expSign * pAnalytical;
pAnalytical = pAnalytical(:);

figure;
hold on;
plot(zPlot * 1e3, abs(pKW), 'LineWidth', 2);
plot(zPlot * 1e3, abs(pAnalytical), '--', 'LineWidth', 2);
xlabel('z, mm');
ylabel('Pressure amplitude, Pa');
title('Water test: k-Wave vs. on-axis analytical (O''Neil)');
legend('k-Wave', 'Analytical');
grid on;
grid minor;
hold off;
end
