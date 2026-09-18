function MaterialMatrixResized = resize_material_matrix_strong(MaterialMatrix, nxPadBegin, nyPadBegin, nzPadBegin, ...
    nxPadEnd, nyPadEnd, nzPadEnd, soundSpeedContactMedium, densityContactMedium, alphaContactMedium)
%RESIZE_MATERIAL_MATRIX_STRONG Resize material fields with minimal temporary copies.
%   This variant is intended for strong memory-saving mode. It allocates the
%   final array for each field once, fills it with the default padding value,
%   and copies only the cropped source block into place.

% Get current dimensions
[nx, ny, nz] = size(MaterialMatrix.c0);

% Calculate new dimensions
nxNew = nx + nxPadBegin + nxPadEnd;
nyNew = ny + nyPadBegin + nyPadEnd;
nzNew = nz + nzPadBegin + nzPadEnd;

% Validate that new dimensions are positive
if nxNew <= 0 || nyNew <= 0 || nzNew <= 0
    error('Resizing would result in non-positive dimensions. nxNew=%d, nyNew=%d, nzNew=%d', ...
        nxNew, nyNew, nzNew);
end

% Preserve only metadata fields without copying the bulky matrices.
MaterialMatrixResized = struct();
for fieldName = ['dx', 'dy', 'dz']
    if isfield(MaterialMatrix, fieldName)
        MaterialMatrixResized.(fieldName) = MaterialMatrix.(fieldName);
    end
end

% List of 3D matrix fields to resize
matrixFields = {'c0', 'rho0', 'alpha'};

% Default values for each field
defaultValues = struct('c0', soundSpeedContactMedium, ...
                      'rho0', densityContactMedium, ...
                      'alpha', alphaContactMedium);

% Source crop limits after trimming
xSrcStart = 1 + max(0, -nxPadBegin);
xSrcEnd = nx - max(0, -nxPadEnd);
ySrcStart = 1 + max(0, -nyPadBegin);
ySrcEnd = ny - max(0, -nyPadEnd);
zSrcStart = 1 + max(0, -nzPadBegin);
zSrcEnd = nz - max(0, -nzPadEnd);

% Destination insertion limits after padding
xDstStart = 1 + max(0, nxPadBegin);
yDstStart = 1 + max(0, nyPadBegin);
zDstStart = 1 + max(0, nzPadBegin);
xDstEnd = xDstStart + (xSrcEnd - xSrcStart);
yDstEnd = yDstStart + (ySrcEnd - ySrcStart);
zDstEnd = zDstStart + (zSrcEnd - zSrcStart);

% Process each 3D matrix field
for i = 1:length(matrixFields)
    fieldName = matrixFields{i};

    if isfield(MaterialMatrix, fieldName)
        sourceMatrix = MaterialMatrix.(fieldName);
        resizedMatrix = defaultValues.(fieldName) * ones(nxNew, nyNew, nzNew, class(sourceMatrix));
        resizedMatrix(xDstStart:xDstEnd, yDstStart:yDstEnd, zDstStart:zDstEnd) = ...
            sourceMatrix(xSrcStart:xSrcEnd, ySrcStart:ySrcEnd, zSrcStart:zSrcEnd);
        MaterialMatrixResized.(fieldName) = resizedMatrix;
    end
end

% Verify final dimensions
if ~isequal(size(MaterialMatrixResized.c0), [nxNew, nyNew, nzNew])
    error('Dimension mismatch after resizing. Expected [%d, %d, %d], got [%d, %d, %d]', ...
        nxNew, nyNew, nzNew, size(MaterialMatrixResized.c0, 1), ...
        size(MaterialMatrixResized.c0, 2), size(MaterialMatrixResized.c0, 3));
end

end
