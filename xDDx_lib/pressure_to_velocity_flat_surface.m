% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

%Power calculation function
function [velocityOut] = pressure_to_velocity_flat_surface(complexPressureAmplitude, stepDim1, stepDim2, frequency, MediumParameters)
waveNumber = 2*pi*frequency/MediumParameters.soundSpeed;

numberOfPointsDim1 = size(complexPressureAmplitude,1);
numberOfPointsDim2 = size(complexPressureAmplitude,2);

[kDim1, ~] = make_k_array(numberOfPointsDim1, stepDim1);
[kDim2, ~] = make_k_array(numberOfPointsDim2, stepDim2);

[Kx, Ky] = meshgrid(kDim2, kDim1);

angularSpectrumPressure = fftshift(fft2(complexPressureAmplitude));

angularSpectrumVelocity = 1/MediumParameters.density/MediumParameters.soundSpeed*sqrt(1-(Kx.^2 + Ky.^2)/waveNumber^2).*angularSpectrumPressure;

angularSpectrumVelocity(Kx.^2 + Ky.^2 >= waveNumber^2) = 0;

velocityOut = ifft2(ifftshift(angularSpectrumVelocity));

end

