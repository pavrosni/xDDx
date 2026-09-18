% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function TransducerSf = read_transducer_sf_from_xls(filePath)

% Use the spreadsheet APIs available in MATLAB R2016b.
[~, sheetNames] = xlsfinfo(filePath);
TransducerSf = [];

% Process scalars
[~, ~, scalarData] = xlsread(filePath, 2, '', 'basic');
for i = 1:size(scalarData, 1)
    fieldName = scalarData{i, 1};
    fieldValue = scalarData{i, 2};

    if ~isvarname(fieldName)
        fieldName = matlab.lang.makeValidName(fieldName);
    end

    if strcmpi(fieldName, 'radiusOfCurvature') && ismissing(fieldValue)
        TransducerSf.(fieldName) = [];
    else
        if isnumeric(fieldValue) && isscalar(fieldValue)
            TransducerSf.(fieldName) = fieldValue;
        else
            error('Invalid value for "%s" in "Scalar Parameters" sheet!', fieldName);
        end
    end
end

% Process x-, y-, z-grids, amplitude, and phase
for sheetIdx = 3:7
    sheetName = sheetNames{sheetIdx};
    matrixData = read_numeric_xls_sheet(filePath, sheetIdx);


    if ~isvarname(sheetName)
        sheetName = matlab.lang.makeValidName(sheetName);
    end

    if sheetIdx < 6
        fieldName = sheetName(1:end-3);
    elseif sheetIdx == 6
        fieldName = sheetName(1:end-7);
    elseif sheetIdx == 7
        fieldName = sheetName(1:end-5);
    end


    if (sheetIdx == 5)  && (isempty(matrixData(:))  || sum(isnan(matrixData(:))))
        if isempty(TransducerSf.radiusOfCurvature)
            matrixData = zeros(size(TransducerSf.xGrid));
        else
            matrixData = TransducerSf.radiusOfCurvature - sqrt(TransducerSf.radiusOfCurvature^2 - TransducerSf.xGrid.^2 - TransducerSf.yGrid.^2);
        end
    end

    if sum(isnan(matrixData(:)))
        error('Invalid element value for "%s" sheet!', fieldName);
    end

    if sheetIdx < 6
        TransducerSf.(fieldName) = matrixData;
    elseif sheetIdx == 6
        TransducerSf.complexVelocityAmplitude = matrixData;
    elseif sheetIdx == 7
        TransducerSf.complexVelocityAmplitude = TransducerSf.complexVelocityAmplitude.*exp(1i*matrixData);
    end
end


end
