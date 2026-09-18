function matrixData = read_numeric_xls_sheet(filePath, sheetIndex)
%READ_NUMERIC_XLS_SHEET Read a numeric worksheet using MATLAB R2016b APIs.
% Preserve invalid cells as NaN for the caller's validation. The numeric
% output of xlsread can silently trim text-only rows and columns.
[~, ~, matrixCells] = xlsread(filePath, sheetIndex, '', 'basic');

% Basic mode can include empty formatted cells beyond the actual data.
% Trim only empty borders, keeping holes within the data for validation.
isEmptyCell = cellfun(@(value) isempty(value) || ...
    (isnumeric(value) && isscalar(value) && isnan(value)), matrixCells);
occupiedRows = find(any(~isEmptyCell, 2));
occupiedColumns = find(any(~isEmptyCell, 1));
if isempty(occupiedRows) || isempty(occupiedColumns)
    matrixData = [];
    return;
end
matrixCells = matrixCells(occupiedRows(1):occupiedRows(end), ...
    occupiedColumns(1):occupiedColumns(end));

% Numeric text is accepted; other text and missing cells remain NaN.
isText = cellfun(@ischar, matrixCells);
matrixCells(isText) = num2cell(str2double(matrixCells(isText)));
isNumeric = cellfun(@(value) isnumeric(value) && isscalar(value), matrixCells);
matrixCells(~isNumeric) = {NaN};
matrixData = cell2mat(matrixCells);
end
