% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [numberOfFrequencySamples] = find_lowest_significant_freq(fieldSpectrum, powerThrsh)

powerAtFrequencyEstimation = squeeze(sum(sum(abs(fieldSpectrum(:,:,2:end)).^2,1),2));

powerAtFrequencyEstimation = powerAtFrequencyEstimation./max(powerAtFrequencyEstimation(:));

numberOfFrequencySamples = find(powerAtFrequencyEstimation > powerThrsh, 1, 'first');

if (numberOfFrequencySamples < 1) && isempty(numberOfFrequencySamples)
  numberOfFrequencySamples = 1;
end

end

