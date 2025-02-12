% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [Geometry, HologramSf, Medium] = read_hologram_sf_from_xls(filePath)

sheetNames = sheetnames(filePath);

% Process Geometry
sheetIdx = 2;
Geometry = [];
scalarData = readcell(filePath, 'Sheet', sheetNames{sheetIdx});
for i = 1:size(scalarData, 1)
    fieldName = scalarData{i, 1};
    fieldValue = scalarData{i, 2};

    if ~isvarname(fieldName)
        fieldName = matlab.lang.makeValidName(fieldName);
    end

    if strcmpi(fieldName, 'radiusOfCurvature') && ismissing(fieldValue)
        Geometry.(fieldName) = [];
    else

        if isnumeric(fieldValue) && isscalar(fieldValue)
            Geometry.(fieldName) = fieldValue;
        else
            error('Invalid value for "%s" in "Geometry" sheet!', fieldName);
        end

    end
end


% Process hologram scalars
sheetIdx = 3;
HologramSf = [];
scalarData = readcell(filePath, 'Sheet', sheetNames{sheetIdx});
for i = 1:size(scalarData, 1)
    fieldName = scalarData{i, 1};
    fieldValue = scalarData{i, 2};

    if ~isvarname(fieldName)
        fieldName = matlab.lang.makeValidName(fieldName);
    end

    if isnumeric(fieldValue) && isscalar(fieldValue)
        HologramSf.(fieldName) = fieldValue;
    else
        error('Invalid value for "%s" in "Scalar Hologram Parameters" sheet!', fieldName);
    end
end


% Process x-, y-, z-grids, amplitude, and phase
for sheetIdx = 4:7
    sheetName = sheetNames{sheetIdx};
    matrixData = readmatrix(filePath, 'Sheet', sheetName);


    if ~isvarname(sheetName)
        sheetName = matlab.lang.makeValidName(sheetName);
    end

    if sheetIdx < 6
        fieldName = sheetName(1:end-3);
    elseif sheetIdx == 6
        fieldName = sheetName(1:end-4);
    elseif sheetIdx == 7
        fieldName = sheetName(1:end-6);
    end

    if sum(isnan(matrixData(:)))
        error('Invalid element value for "%s" Sheet!', fieldName);
    end

    if sheetIdx < 6
        HologramSf.(fieldName) = matrixData;
    elseif sheetIdx == 6
        HologramSf.complexPressureAmplitude = matrixData;
    elseif sheetIdx == 7
        HologramSf.complexPressureAmplitude = HologramSf.complexPressureAmplitude.*exp(1i*matrixData);
    end
end

% Process hologram scalars
sheetIdx = 8;
Medium = [];
scalarData = readcell(filePath, 'Sheet', sheetNames{sheetIdx});
for i = 1:size(scalarData, 1)
    fieldName = scalarData{i, 1};
    fieldValue = scalarData{i, 2};

    if ~isvarname(fieldName)
        fieldName = matlab.lang.makeValidName(fieldName);
    end

    if isnumeric(fieldValue) && isscalar(fieldValue)
        Medium.(fieldName) = fieldValue;
    else
        error('Invalid value for "%s" in "Medium" sheet!', fieldName);
    end
end


end
