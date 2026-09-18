function TransducerSf = load_xDDx_transducer(dataPath)
%LOAD_XDDX_TRANSDUCER Load and validate TransducerSf from a MAT file.
%
%   TransducerSf = load_xDDx_transducer(dataPath)
%
%   The MAT file must contain a variable named TransducerSf (struct) matching
%   the single-frequency transducer format used by the xDDx Rayleigh toolbox
%   (see transducer_simulation_sf.m header).
%
%   Required fields:
%     expSign, frequency, xGrid, yGrid, dx, dy, complexVelocityAmplitude
%
%   Optional fields:
%     radiusOfCurvature  - [] or absent for flat; positive scalar for spherical
%     zGrid              - if omitted, set to zeros (flat) or sphere surface
%                          z = R - sqrt(R^2 - x^2 - y^2) when R is given
%     type               - if present, must be 'custom'. The returned
%                          TransducerSf.type is always set to 'custom'.

    if nargin < 1 || isempty(dataPath)
        error('dataPath must be a non-empty char vector or string.');
    end
    if isa(dataPath, 'string')
        dataPath = char(dataPath);
    end
    if ~ischar(dataPath)
        error('dataPath must be a char vector or string.');
    end

    if exist(dataPath, 'file') ~= 2
        error('File not found: %s', dataPath);
    end

    S = load(dataPath);
    if ~isfield(S, 'TransducerSf')
        vars = strjoin(fieldnames(S), ', ');
        error('MAT file must contain variable ''TransducerSf''. Found: %s', vars);
    end

    TransducerSf = validate_TransducerSf_structure(S.TransducerSf);
end

function TransducerSf = validate_TransducerSf_structure(T)
    if ~isstruct(T) || numel(T) ~= 1
        error('TransducerSf must be a scalar struct.');
    end

    if isfield(T, 'type')
        inputType = T.type;
        if isa(inputType, 'string')
            inputType = char(inputType);
        end
        if ~ischar(inputType) || ~strcmp(inputType, 'custom')
            error('TransducerSf.type from load_xDDx_transducer must be ''custom'' if present. Delete the field or set it to ''custom''.');
        end
    end

    required = {'expSign', 'frequency', 'xGrid', 'yGrid', 'dx', 'dy', 'complexVelocityAmplitude'};
    for k = 1:numel(required)
        if ~isfield(T, required{k}) || isempty(T.(required{k}))
            error('TransducerSf missing or empty required field: %s', required{k});
        end
    end

    expSign = T.expSign;
    if ~isnumeric(expSign) || ~isscalar(expSign) || ~any(expSign == [-1, 1])
        error('TransducerSf.expSign must be scalar +1 or -1.');
    end

    frequency = T.frequency;
    if ~isnumeric(frequency) || ~isscalar(frequency) || frequency <= 0
        error('TransducerSf.frequency must be a positive scalar (Hz).');
    end

    dx = T.dx;
    dy = T.dy;
    if ~isnumeric(dx) || ~isscalar(dx) || dx <= 0
        error('TransducerSf.dx must be a positive scalar.');
    end
    if ~isnumeric(dy) || ~isscalar(dy) || dy <= 0
        error('TransducerSf.dy must be a positive scalar.');
    end

    xGrid = T.xGrid;
    yGrid = T.yGrid;
    v = T.complexVelocityAmplitude;
    if ~isnumeric(xGrid) || ~isnumeric(yGrid) || ~isnumeric(v)
        error('TransducerSf.xGrid, yGrid and complexVelocityAmplitude must be numeric.');
    end
    if ~isequal(size(xGrid), size(yGrid), size(v))
        error('TransducerSf: xGrid, yGrid and complexVelocityAmplitude must have the same size.');
    end

    if isfield(T, 'radiusOfCurvature') && ~isempty(T.radiusOfCurvature)
        R = T.radiusOfCurvature;
        if ~isnumeric(R) || ~isscalar(R) || R <= 0
            error('TransducerSf.radiusOfCurvature must be a positive scalar or [].');
        end
    else
        R = [];
    end

    if isfield(T, 'zGrid') && ~isempty(T.zGrid)
        zGrid = T.zGrid;
        if ~isnumeric(zGrid)
            error('TransducerSf.zGrid must be numeric.');
        end
        if ~isequal(size(zGrid), size(xGrid))
            error('TransducerSf.zGrid must match the size of xGrid and yGrid.');
        end
    else
        if isempty(R)
            zGrid = zeros(size(xGrid));
        else
            rsq = xGrid.^2 + yGrid.^2;
            zGrid = R - sqrt(max(0, R^2 - rsq));
        end
    end

    TransducerSf = [];
    TransducerSf.type = 'custom';
    TransducerSf.expSign = double(expSign);
    TransducerSf.frequency = double(frequency);
    TransducerSf.xGrid = xGrid;
    TransducerSf.yGrid = yGrid;
    TransducerSf.zGrid = zGrid;
    TransducerSf.dx = double(dx);
    TransducerSf.dy = double(dy);
    TransducerSf.complexVelocityAmplitude = v;
    if isempty(R)
        TransducerSf.radiusOfCurvature = [];
    else
        TransducerSf.radiusOfCurvature = double(R);
    end
end
