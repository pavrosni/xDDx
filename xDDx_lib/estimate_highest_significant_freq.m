% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [numberOfFrequencySamples] = estimate_highest_significant_freq(frequencyStep, pSpectrum, powerThrsh, maxSamples, maxFrequency)

powerAtFrequencyEstimation = zeros(1,(size(pSpectrum,3)-1));
for iFreq = 1:length(powerAtFrequencyEstimation)
    powerAtFrequencyEstimation(iFreq) = sum(sum(abs(pSpectrum(:,:,iFreq+1)).^2,1),2);
end

powerAtFrequencyEstimation = powerAtFrequencyEstimation./max(powerAtFrequencyEstimation(:));

numberOfFrequencySamples = find(powerAtFrequencyEstimation > powerThrsh, 1, 'last')+1;

if (numberOfFrequencySamples > maxSamples) || (frequencyStep*numberOfFrequencySamples > maxFrequency)
    warning('The number of frequency samples that fulfills the threshold criteria (numberOfFrequencySamples) is greater than the maximum number of frequency samples (maxSamples). To avoid performance issues, numberOfFrequencySamples is set equal to maxSamples.');
    numberOfFrequencySamples = min(maxSamples, fix(maxFrequency/frequencyStep)+1);
end
end

