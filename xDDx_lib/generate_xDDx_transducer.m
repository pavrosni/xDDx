function TransducerSf = generate_xDDx_transducer(varargin)
%GENERATE_XDDX_TRANSDUCER Generate flat or spherical TransducerSf.
%   Backward-compatible positional usage:
%     generate_xDDx_transducer(aperture, radiusOfCurvature, frequency, dx, dy, initialVelocity)
%
%   Single-element name-value usage:
%     generate_xDDx_transducer('aperture', 50e-3, ...
%         'frequency', 1e6, 'initialVelocity', 1, ...)
%
%   Multi-element array name-value usage:
%     generate_xDDx_transducer('aperture', 50e-3, 'radiusOfCurvature', 50e-3, ...
%         'frequency', 1e6, 'elementAperture', 5e-3, ...
%         'xCenters', xCenters, 'yCenters', yCenters, ...
%         'initialVelocityAmplitudes', amplitudes, 'initialVelocityPhases', phases, ...)
%
%   Required name-value parameters for all transducers:
%     'aperture'            (required) aperture diameter in m
%     'frequency'           (required) frequency in Hz
%
%   Required name-value parameter for single-element transducers:
%     'initialVelocity'     (required) velocity amplitude in m/s
%
%   Optional name-value parameters:
%     'radiusOfCurvature'   spherical curvature radius in m (default [] for flat)
%     'sourceStepX'         x-step in m (default 3 points per water wavelength)
%     'sourceStepY'         y-step in m (default 3 points per water wavelength)
%     'elementAperture'     circular array-element aperture in m
%     'xCenters'            x-coordinates of array-element centers in m
%     'yCenters'            y-coordinates of array-element centers in m
%     'initialVelocityAmplitudes' array-element velocity amplitudes in m/s
%     'initialVelocityPhases'     array-element velocity phases in rad
%     'expSign'             +1 or -1 (default 1)
%
%   The five array-element parameters are optional as a group. If any one is
%   provided, all five must be provided.
%
%   The default optional parameters are:
%     'radiusOfCurvature' = []
%     'sourceStepX' = waterSoundSpeed/frequency/3
%     'sourceStepY' = waterSoundSpeed/frequency/3
%     'expSign' = 1

    if nargin == 6 && isnumeric(varargin{1})
        params = [];
        params.aperture = varargin{1};
        params.radiusOfCurvature = varargin{2};
        params.frequency = varargin{3};
        params.sourceStepX = varargin{4};
        params.sourceStepY = varargin{5};
        params.initialVelocity = varargin{6};
        params.elementAperture = [];
        params.xCenters = [];
        params.yCenters = [];
        params.initialVelocityAmplitudes = [];
        params.initialVelocityPhases = [];
        params.expSign = 1;
        isArrayTransducer = false;
    else
        p = inputParser;
        p.addParameter('aperture', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
        p.addParameter('frequency', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
        p.addParameter('sourceStepX', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
        p.addParameter('sourceStepY', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
        p.addParameter('initialVelocity', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
        p.addParameter('radiusOfCurvature', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
        p.addParameter('elementAperture', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
        p.addParameter('xCenters', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
        p.addParameter('yCenters', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
        p.addParameter('initialVelocityAmplitudes', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
        p.addParameter('initialVelocityPhases', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
        p.addParameter('expSign', 1, @(x) isnumeric(x) && isscalar(x) && any(x == [-1, 1]));
        p.parse(varargin{:});
        params = p.Results;

        arrayParameterNames = {'elementAperture', 'xCenters', 'yCenters', ...
            'initialVelocityAmplitudes', 'initialVelocityPhases'};
        isArrayParameterProvided = ~ismember(arrayParameterNames, p.UsingDefaults);
        isArrayTransducer = any(isArrayParameterProvided);
        if isArrayTransducer && ~all(isArrayParameterProvided)
            missingArrayParameters = arrayParameterNames(~isArrayParameterProvided);
            error('Array transducer parameters must be provided together. Missing: %s.', ...
                strjoin(missingArrayParameters, ', '));
        end

        requiredParameters = {
            'aperture', 'aperture is not set, a value is needed (e.g. ''aperture'', 50e-3)';
            'frequency', 'frequency is not set, a value is needed (e.g. ''frequency'', 1e6)';
            'initialVelocity', 'initialVelocity is not set, a value is needed (e.g. ''initialVelocity'', 1/1000/1500)';
        };
        if isArrayTransducer
            requiredParameters(strcmp(requiredParameters(:, 1), 'initialVelocity'), :) = [];
        end
        for iParameter = 1:size(requiredParameters, 1)
            if any(strcmp(p.UsingDefaults, requiredParameters{iParameter, 1}))
                error(requiredParameters{iParameter, 2});
            end
        end

        waterSoundSpeed = 1500; % m/s
        defaultSourceStep = waterSoundSpeed / params.frequency / 3;
        if isempty(params.sourceStepX)
            params.sourceStepX = defaultSourceStep;
        end
        if isempty(params.sourceStepY)
            params.sourceStepY = defaultSourceStep;
        end
    end

    if isempty(params.aperture)
        error('aperture must be a positive scalar.');
    end
    if isempty(params.frequency)
        error('frequency must be a positive scalar.');
    end
    if ~isArrayTransducer && isempty(params.initialVelocity)
        error('initialVelocity must be a scalar.');
    end

    if isArrayTransducer
        if isempty(params.radiusOfCurvature)
            error('radiusOfCurvature is needed for array transducer generation (e.g. ''radiusOfCurvature'', 50e-3).');
        end
        nElements = numel(params.xCenters);
        if nElements == 0
            error('Array transducer parameters must define at least one element.');
        end
        if numel(params.yCenters) ~= nElements || ...
                numel(params.initialVelocityAmplitudes) ~= nElements || ...
                numel(params.initialVelocityPhases) ~= nElements
            error('xCenters, yCenters, initialVelocityAmplitudes, and initialVelocityPhases must have the same number of elements.');
        end
        if any(params.xCenters(:).^2 + params.yCenters(:).^2 > (params.aperture / 2)^2)
            error('Array element centers must lie inside the array aperture: xCenters.^2 + yCenters.^2 must be <= (aperture/2)^2.');
        end
        if params.elementAperture >= 2 * params.radiusOfCurvature
            error('Array element aperture is impossible: elementAperture must be < 2*radiusOfCurvature.');
        end
        if any(params.xCenters(:).^2 + params.yCenters(:).^2 >= params.radiusOfCurvature^2)
            error('Array element centers must lie inside the spherical surface projection: xCenters.^2 + yCenters.^2 must be < radiusOfCurvature^2.');
        end
    end

    isSpherical = ~isempty(params.radiusOfCurvature);
    if isSpherical && (params.aperture >= 2 * params.radiusOfCurvature)
        error('Spherical transducer is impossible: aperture must be < 2*radiusOfCurvature.');
    end

    nxSource = round(params.aperture / params.sourceStepX) + 2;
    nySource = round(params.aperture / params.sourceStepY) + 2;
    [xSource, ySource] = build_flat_grid_centered(nxSource, params.sourceStepX, nySource, params.sourceStepY);

    if isSpherical
        zSource = params.radiusOfCurvature - sqrt(params.radiusOfCurvature^2 - xSource.^2 - ySource.^2);
        radiusOfCurvature = params.radiusOfCurvature;
    else
        zSource = zeros(size(xSource));
        radiusOfCurvature = [];
    end

    if isArrayTransducer
        complexVelocityAmplitude = zeros(size(xSource));
        elementSectionRadius = params.elementAperture / 2;
        elementCosHalfAngle = sqrt(params.radiusOfCurvature^2 - elementSectionRadius^2) / params.radiusOfCurvature;
        xCenters = params.xCenters(:);
        yCenters = params.yCenters(:);
        initialVelocityAmplitudes = params.initialVelocityAmplitudes(:);
        initialVelocityPhases = params.initialVelocityPhases(:);
        for iElement = 1:nElements
            zCenter = params.radiusOfCurvature - sqrt(params.radiusOfCurvature^2 ...
                - xCenters(iElement)^2 - yCenters(iElement)^2);
            sphereDotProduct = xSource .* xCenters(iElement) ...
                + ySource .* yCenters(iElement) ...
                + (zSource - params.radiusOfCurvature) .* (zCenter - params.radiusOfCurvature);
            elementMask = sphereDotProduct >= params.radiusOfCurvature^2 * elementCosHalfAngle;
            complexVelocityAmplitude(elementMask) = initialVelocityAmplitudes(iElement) ...
                * exp(1i * initialVelocityPhases(iElement));
        end
    else
        complexVelocityAmplitude = params.initialVelocity * ones(size(xSource));
    end
    complexVelocityAmplitude(xSource.^2 + ySource.^2 > (params.aperture / 2)^2) = 0;

    TransducerSf = [];
    if isArrayTransducer
        TransducerSf.type = 'standard_array';
    else
        TransducerSf.type = 'standard_single_element';
    end
    TransducerSf.expSign = params.expSign;
    TransducerSf.frequency = params.frequency;
    TransducerSf.radiusOfCurvature = radiusOfCurvature;
    TransducerSf.xGrid = xSource;
    TransducerSf.yGrid = ySource;
    TransducerSf.zGrid = zSource;
    TransducerSf.dx = params.sourceStepX;
    TransducerSf.dy = params.sourceStepY;
    TransducerSf.complexVelocityAmplitude = complexVelocityAmplitude;
end
