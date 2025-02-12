% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [numberOfFrequencySamples] = find_highest_significant_freq(frequencyStep, fieldSpectrum, powerThrsh, maxSamples, maxFrequency)

powerAtFrequencyEstimation = squeeze(sum(sum(abs(fieldSpectrum(:,:,2:end)).^2,1),2));

powerAtFrequencyEstimation = powerAtFrequencyEstimation./max(powerAtFrequencyEstimation(:));

numberOfFrequencySamples = find(powerAtFrequencyEstimation > powerThrsh, 1, 'last')+1;

if (numberOfFrequencySamples > maxSamples) || (frequencyStep*numberOfFrequencySamples > maxFrequency)
    warning('The number of frequency samples that fulfills the threshold criteria (maxFrequencySample) is greater than the maximum number of frequency samples (maxNumFrequencySamples). To avoid performance issues, numberOfFrequencySamples is set equal to maxSamples.');
    numberOfFrequencySamples = min(maxSamples, fix(maxFrequency/frequencyStep)+1);
end
end

