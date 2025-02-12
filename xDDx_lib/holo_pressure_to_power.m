% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [powerOut] = holo_pressure_to_power(complexPressureAmplitude, stepDim1, stepDim2, frequency, Medium)
waveNumber = 2*pi*frequency/Medium.soundSpeed;

numberOfPointsDim1 = size(complexPressureAmplitude,1);
numberOfPointsDim2 = size(complexPressureAmplitude,2);

[kDim1, kStepDim1] = make_k_array(numberOfPointsDim1, stepDim1);
[kDim2, kStepDim2] = make_k_array(numberOfPointsDim2, stepDim2);

[Kx, Ky] = meshgrid(kDim2, kDim1);

angularSpectrum = fftshift(fft2(complexPressureAmplitude)*stepDim1*stepDim2);

intensity = abs(angularSpectrum).^2.*sqrt(1 - (Kx.^2 + Ky.^2)/waveNumber^2)*kStepDim1*kStepDim2;
intensityCircle = intensity(Kx.^2 + Ky.^2 < waveNumber^2);

powerOut = 1/(8*pi^2*Medium.density*Medium.soundSpeed)*sum(sum(intensityCircle));

end

