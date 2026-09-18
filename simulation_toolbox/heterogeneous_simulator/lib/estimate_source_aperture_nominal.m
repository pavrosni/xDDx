function [apertureNominalX, apertureNominalY] = estimate_source_aperture_nominal(TransducerSf, thresholdFraction)
%ESTIMATE_SOURCE_APERTURE_NOMINAL Estimate active source aperture in x and y.
%   [apertureNominalX, apertureNominalY] = estimate_source_aperture_nominal(TransducerSf)
%   estimates the xy projection of the active source aperture using points
%   with velocity magnitude greater than 1% of the maximum magnitude.

    if nargin < 2 || isempty(thresholdFraction)
        thresholdFraction = 0.01;
    end
    if ~isnumeric(thresholdFraction) || ~isscalar(thresholdFraction) || ...
            thresholdFraction < 0 || thresholdFraction >= 1
        error('thresholdFraction must be a scalar in the range [0, 1).');
    end

    requiredFields = {'xGrid', 'yGrid', 'complexVelocityAmplitude'};
    for iField = 1:numel(requiredFields)
        if ~isfield(TransducerSf, requiredFields{iField})
            error('TransducerSf.%s is required to estimate nominal aperture.', requiredFields{iField});
        end
    end

    xGrid = TransducerSf.xGrid;
    yGrid = TransducerSf.yGrid;
    velocityMagnitude = abs(TransducerSf.complexVelocityAmplitude);
    if ~isequal(size(xGrid), size(yGrid), size(velocityMagnitude))
        error('TransducerSf.xGrid, yGrid, and complexVelocityAmplitude must have the same size.');
    end

    maxVelocityMagnitude = max(velocityMagnitude(:));
    if isempty(maxVelocityMagnitude) || maxVelocityMagnitude <= 0
        error('Cannot estimate nominal aperture because TransducerSf.complexVelocityAmplitude is zero everywhere.');
    end

    activeMask = velocityMagnitude > thresholdFraction * maxVelocityMagnitude;
    if ~any(activeMask(:))
        error('Cannot estimate nominal aperture because no source points exceed the active threshold.');
    end

    activeX = xGrid(activeMask);
    activeY = yGrid(activeMask);
    apertureNominalX = max(activeX(:)) - min(activeX(:));
    apertureNominalY = max(activeY(:)) - min(activeY(:));

    if apertureNominalX <= 0 || apertureNominalY <= 0
        error('Estimated nominal aperture must be positive in both x and y directions.');
    end
end
