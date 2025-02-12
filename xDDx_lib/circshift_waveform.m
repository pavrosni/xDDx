% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [vWaveformsSurfOut, iBegin, iEnd] = circshift_waveform(vWaveformsSurf,levelSignalToSave,showRingingDown)
    
    maxSgnlTime= squeeze(max(max(abs(vWaveformsSurf),[],1),[],2));
    [~,idxtMaxSgnl] = max(maxSgnlTime);
    iMiddle = fix(size(vWaveformsSurf,3)/2);
    vWaveformsSurfCentered = circshift(vWaveformsSurf,(iMiddle-idxtMaxSgnl-1),3);
    maxSgnlTimeCentered = circshift(maxSgnlTime,(iMiddle-idxtMaxSgnl-1));

    maxSgnlTimeFirstHalfFlipped = flip(maxSgnlTimeCentered(1:iMiddle));

    iBeginHalfFlipped = find(maxSgnlTimeFirstHalfFlipped < levelSignalToSave * max(maxSgnlTime), 1, 'first');

    if isempty(iBeginHalfFlipped)
      [~,iBeginHalfFlipped] = min(maxSgnlTimeFirstHalfFlipped);
    end

    iBegin = max(iMiddle-iBeginHalfFlipped, 1);


    if ~showRingingDown
        iEnd = find(maxSgnlTimeCentered > levelSignalToSave * max(maxSgnlTimeCentered), 1, 'last');

        if isempty(iEnd)
            iEnd = length(maxSgnlTimeCentered);
        end
    else
        iEnd = length(maxSgnlTimeCentered);
    end

    vWaveformsSurfOut = vWaveformsSurfCentered;
end

