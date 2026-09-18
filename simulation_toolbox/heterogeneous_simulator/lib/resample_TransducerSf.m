function TransducerSf = resample_TransducerSf(TransducerSf, ...
                                        xFlatBoundary, yFlatBoundary, zFlatBoundary, ...
                                        dxFlatBoundary, dyFlatBoundary)
%RESAMPLE_TRANSDUCERSF Resample TransducerSf by nearest-neighbor interpolation.
%
%   TransducerSf = resample_TransducerSf(TransducerSf, xFlatBoundary,
%   yFlatBoundary, zFlatBoundary, dxFlatBoundary, dyFlatBoundary)
%
%   The function maps TransducerSf.complexVelocityAmplitude from the
%   current TransducerSf.xGrid / yGrid to xFlatBoundary / yFlatBoundary
%   using nearest-neighbor interpolation.
%
%   If target grid points already coincide with the transducer grid,
%   interpolation is skipped and no warning is shown.
%
%   NOTE:
%   If interpolation is used, the function warns that users should ideally
%   prepare TransducerSf directly on the MaterialMatrix grid to avoid
%   precision issues.

    requiredFields = {'expSign','frequency','xGrid','yGrid','zGrid','dx','dy','complexVelocityAmplitude'};
    for iField = 1:numel(requiredFields)
        if ~isfield(TransducerSf, requiredFields{iField})
            error('TransducerSf.%s is required.', requiredFields{iField});
        end
    end

    if isfield(TransducerSf, 'radiusOfCurvature') && ~isempty(TransducerSf.radiusOfCurvature)
        error('resample_TransducerSf supports only flat transducers (empty/absent radiusOfCurvature).');
    end
    if ~isequal(size(xFlatBoundary), size(yFlatBoundary), size(zFlatBoundary))
        error('xFlatBoundary, yFlatBoundary and zFlatBoundary must have identical sizes.');
    end

    xGrid0 = TransducerSf.xGrid;
    yGrid0 = TransducerSf.yGrid;
    vGrid0 = TransducerSf.complexVelocityAmplitude;
    zGrid0 = TransducerSf.zGrid;

    if ~isequal(size(xGrid0), size(yGrid0), size(zGrid0), size(vGrid0))
        error('xGrid, yGrid, zGrid and complexVelocityAmplitude must have identical sizes.');
    end

    gridsCoincide = isequal(size(xGrid0), size(xFlatBoundary)) && ...
        (max(abs(xGrid0(:) - xFlatBoundary(:))) < 10 * eps('single')) && ...
        (max(abs(yGrid0(:) - yFlatBoundary(:))) < 10 * eps('single'));

    if gridsCoincide
        vResampled = vGrid0;
    else
        warning(['TransducerSf grid was interpolated onto the simulation boundary grid. ' ...
                 'It is recommended to prepare TransducerSf directly on the same grid ' ...
                 'as MaterialMatrix to avoid possible precision issues.']);

        activeMask = abs(vGrid0) > eps('single');
        if ~any(activeMask(:))
            error(['All transducer grid points have negligible complexVelocityAmplitude ', ...
                '(|v| <= eps(''single'')).']);
        end
        Fv = scatteredInterpolant(xGrid0(activeMask), yGrid0(activeMask), vGrid0(activeMask), 'nearest', 'none');
        vResampled = Fv(xFlatBoundary, yFlatBoundary);
        vResampled(~isfinite(vResampled)) = 0;
    end

    TransducerSf.xGrid = xFlatBoundary;
    TransducerSf.yGrid = yFlatBoundary;
    TransducerSf.zGrid = zFlatBoundary;
    TransducerSf.dx = dxFlatBoundary;
    TransducerSf.dy = dyFlatBoundary;
    TransducerSf.complexVelocityAmplitude = vResampled;
    if isfield(TransducerSf, 'radiusOfCurvature')
        TransducerSf.radiusOfCurvature = [];
    end
end
