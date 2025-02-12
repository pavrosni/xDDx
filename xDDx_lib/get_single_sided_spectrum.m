% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [fArray, fieldSpectrum] = get_single_sided_spectrum(time, waveforms)

%Perform fft to ex tract spectral components for signals at all hologram points 
timeStep = time(2) - time(1);
numberOfTimePoints = size(waveforms,3);
frequencyStep = (timeStep*numberOfTimePoints)^-1;

fieldSpectrum = fft(waveforms,[],3);

fieldSpectrum = fieldSpectrum/numberOfTimePoints;
if mod(numberOfTimePoints,2) == 0    
    fieldSpectrum = fieldSpectrum(:,:,1:(numberOfTimePoints/2+1));
    fieldSpectrum(:,:,2:(end-1)) = 2*fieldSpectrum(:,:,2:(end-1));    
else    
    fieldSpectrum = fieldSpectrum(:,:,1:((numberOfTimePoints+1)/2));
    fieldSpectrum(:,:,2:end) = 2*fieldSpectrum(:,:,2:end);    
end

%Extract a spectral component closest to the frequencyNominal 
fArray = (0:(size(fieldSpectrum,3)-1))*frequencyStep;

end

