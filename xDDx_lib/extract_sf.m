% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [HologramSf, frequency] = extract_sf(HologramTr,expSign,frequencyNominal)

%Perform fft to extract spectral components for signals at all hologram points
timeStep = HologramTr.time(2) - HologramTr.time(1);
numberOfTimePoints = size(HologramTr.pressureWaveforms,3);
frequencyStep = (timeStep*numberOfTimePoints)^-1;

pSpectrum = fft(HologramTr.pressureWaveforms,[],3);

pSpectrum = pSpectrum/numberOfTimePoints;
if mod(numberOfTimePoints,2) == 0
    pSpectrum = pSpectrum(:,:,1:(numberOfTimePoints/2+1));
    pSpectrum(:,:,2:(end-1)) = 2*pSpectrum(:,:,2:(end-1));
else
    pSpectrum = pSpectrum(:,:,1:((numberOfTimePoints+1)/2));
    pSpectrum(:,:,2:end) = 2*pSpectrum(:,:,2:end);
end

if expSign == -1
    pSpectrum = conj(pSpectrum);
end

%Extract a spectral component closest to the frequencyNominal
fArray = (0:(size(pSpectrum,3)-1))*frequencyStep;
[~, ifNominal] = min(abs(fArray - frequencyNominal));
frequency = fArray(ifNominal);

HologramSf = [];
HologramSf.expSign = expSign;
HologramSf.frequency = frequency;
HologramSf.xGrid = HologramTr.xGrid;
HologramSf.yGrid = HologramTr.yGrid;
HologramSf.zPosition = HologramTr.zPosition;
HologramSf.dx = HologramTr.dx;
HologramSf.dy = HologramTr.dy;
HologramSf.complexPressureAmplitude = squeeze((pSpectrum(:,:, ifNominal)));

end

