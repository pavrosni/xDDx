% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [saveSurfSpectrumUpdated, saveSurfSignalUpdated] = check_octave_based_video_issues(isOctave, saveSurfSpectrum, saveSurfSignal)
saveSurfSpectrumUpdated = saveSurfSpectrum;
saveSurfSignalUpdated = saveSurfSignal;
if isOctave && (saveSurfSpectrum || saveSurfSignal)
    warning('Octave does not support ''VideoWriter,'' so the video cannot be recorded! The results will be presented in a frame-by-frame format only.');

    saveSurfSpectrumUpdated = false;
    saveSurfSignalUpdated = false;
end
end

