function TransducerSf = resample_TransducerSf_via_rayleigh(TransducerSf, soundSpeed, density, simulationDevice, ...
                                        xFlatBoundary, yFlatBoundary, zFlatBoundary, ...
                                        dxFlatBoundary, dyFlatBoundary, ...
                                        shiftHoloDistInWl, ServiceParameters)
%RESAMPLE_TRANSDUCERSF_VIA_RAYLEIGH Resample a flat transducer using Rayleigh projections.
%
%   TransducerSf = resample_TransducerSf_via_rayleigh(TransducerSf, soundSpeed, density,
%   simulationDevice, xFlatBoundary, yFlatBoundary, zFlatBoundary,
%   dxFlatBoundary, dyFlatBoundary, shiftHoloDistInWl, ServiceParameters)
%
%   The function is intended for a flat transducer representation
%   (TransducerSf.radiusOfCurvature is empty or absent). The input velocity
%   on the transducer plane is:
%     1) Forward-projected (regime 4: V plane -> P points) to a parallel
%        plane shifted by shiftHoloDist.
%     2) Back-projected (regime 1: P plane -> V plane) to the original
%        transducer plane on a new Cartesian sampling grid.
%
%   The output/resampling grid is taken directly from the provided flat
%   boundary grids xFlatBoundary, yFlatBoundary, zFlatBoundary.

    if nargin < 11 || isempty(ServiceParameters)
        ServiceParameters = [];
        ServiceParameters.threadsPerBlockGPU = 128;
    end

    if nargin < 10 || isempty(shiftHoloDistInWl)
        shiftHoloDistInWl = 10;
    end

    requiredFields = {'expSign','frequency','xGrid','yGrid','zGrid','dx','dy','complexVelocityAmplitude'};
    for iField = 1:numel(requiredFields)
        if ~isfield(TransducerSf, requiredFields{iField})
            error('TransducerSf.%s is required.', requiredFields{iField});
        end
    end

    if isfield(TransducerSf, 'radiusOfCurvature') && ~isempty(TransducerSf.radiusOfCurvature)
        error('resample_TransducerSf_via_rayleigh supports only flat transducers (empty/absent radiusOfCurvature).');
    end
    if ~isequal(size(xFlatBoundary), size(yFlatBoundary), size(zFlatBoundary))
        error('xFlatBoundary, yFlatBoundary and zFlatBoundary must have identical sizes.');
    end

    Medium = [];
    Medium.soundSpeed = soundSpeed;
    Medium.density = density;

    expSign = TransducerSf.expSign;
    frequency = TransducerSf.frequency;
    lambda = Medium.soundSpeed / frequency;
    shiftHoloDist = shiftHoloDistInWl * lambda;

    xGrid0 = TransducerSf.xGrid;
    yGrid0 = TransducerSf.yGrid;
    zGrid0 = TransducerSf.zGrid;
    vGrid0 = TransducerSf.complexVelocityAmplitude;

    if ~isequal(size(xGrid0), size(yGrid0), size(zGrid0), size(vGrid0))
        error('xGrid, yGrid, zGrid and complexVelocityAmplitude must have identical sizes.');
    end

    zMean = mean(zGrid0(:));
    if max(abs(zGrid0(:) - zMean)) > 10 * eps(max(1, abs(zMean)))
        error('Flat transducer expected: TransducerSf.zGrid must lie on one plane.');
    end

    zBoundaryMean = mean(zFlatBoundary(:));
    if max(abs(zFlatBoundary(:) - zBoundaryMean)) > 10 * eps(max(1, abs(zBoundaryMean)))
        error('Flat boundary expected: zFlatBoundary must lie on one plane.');
    end

    % Use the provided flat boundary grid directly.
    xResampled = xFlatBoundary;
    yResampled = yFlatBoundary;
    zResampled = zFlatBoundary;

    xShift = xResampled;
    yShift = yResampled;
    zShift = zResampled + shiftHoloDist;

    SourceParameters = [];
    withinTheActiveSurface = (abs(vGrid0) > eps);
    SourceParameters.xGrid = xGrid0(withinTheActiveSurface);
    SourceParameters.yGrid = yGrid0(withinTheActiveSurface);
    SourceParameters.zGrid = zGrid0(withinTheActiveSurface);
    SourceParameters.dx = TransducerSf.dx;
    SourceParameters.dy = TransducerSf.dy;
    SourceParameters.input = vGrid0(withinTheActiveSurface);

    FieldParameters = [];
    FieldParameters.xGrid = xShift;
    FieldParameters.yGrid = yShift;
    FieldParameters.zGrid = zShift;

    regime = 4; % Forward-projection: V on a plane --> P at points
    isTransient = false;
    pShift = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, ...
        SourceParameters, FieldParameters, Medium, ServiceParameters);

    SourceParameters = [];
    SourceParameters.xGrid = xResampled;
    SourceParameters.yGrid = yResampled;
    SourceParameters.zGrid = zResampled;

    FieldParameters = [];
    FieldParameters.xGrid = xShift;
    FieldParameters.yGrid = yShift;
    FieldParameters.zGrid = zShift;
    FieldParameters.dx = dxFlatBoundary;
    FieldParameters.dy = dyFlatBoundary;
    FieldParameters.input = pShift;

    regime = 1; % Back-projection: P on a plane --> V on a plane
    vResampled = rayleigh_simulator(expSign, frequency, regime, simulationDevice, isTransient, ...
        SourceParameters, FieldParameters, Medium, ServiceParameters);

    TransducerSf.xGrid = xResampled;
    TransducerSf.yGrid = yResampled;
    TransducerSf.zGrid = zResampled;
    TransducerSf.dx = dxFlatBoundary;
    TransducerSf.dy = dyFlatBoundary;
    TransducerSf.complexVelocityAmplitude = vResampled;
    if isfield(TransducerSf, 'radiusOfCurvature')
        TransducerSf.radiusOfCurvature = [];
    end
end
