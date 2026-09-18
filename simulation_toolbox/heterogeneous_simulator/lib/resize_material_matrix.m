function MaterialMatrixResized = resize_material_matrix(MaterialMatrix, nxPadBegin, nyPadBegin, nzPadBegin, ...
    nxPadEnd, nyPadEnd, nzPadEnd, soundSpeedContactMedium, densityContactMedium, alphaContactMedium)
%RESIZE_MATERIAL_MATRIX Resize MaterialMatrix structure based on padding values.
%   MaterialMatrixResized = resize_material_matrix(MaterialMatrix, nxPadBegin, nyPadBegin, nzPadBegin, ...
%       nxPadEnd, nyPadEnd, nzPadEnd, soundSpeedContactMedium, densityContactMedium, alphaContactMedium)
%
%   This function resizes all 3D matrix fields in MaterialMatrix (c0, rho0, alpha)
%   based on the padding values:
%   - Positive values: Add default-value elements at the beginning/end of the dimension
%   - Negative values: Remove that number of elements from the beginning/end of the dimension
%
%   Inputs:
%       MaterialMatrix - Structure containing 3D matrices (c0, rho0, alpha) and scalars (dx, dy, dz)
%       nxPadBegin - Number of elements to add (positive) or remove (negative) at beginning of x dimension
%       nyPadBegin - Number of elements to add (positive) or remove (negative) at beginning of y dimension
%       nzPadBegin - Number of elements to add (positive) or remove (negative) at beginning of z dimension
%       nxPadEnd - Number of elements to add (positive) or remove (negative) at end of x dimension
%       nyPadEnd - Number of elements to add (positive) or remove (negative) at end of y dimension
%       nzPadEnd - Number of elements to add (positive) or remove (negative) at end of z dimension
%       soundSpeedContactMedium - Default value for c0 when padding (default: 1500)
%       densityContactMedium - Default value for rho0 when padding (default: 1000)
%       alphaContactMedium - Default value for alpha when padding (default: 0)
%
%   Output:
%       MaterialMatrixResized - Resized MaterialMatrix structure

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

% Initialize output structure
MaterialMatrixResized = MaterialMatrix;

% List of 3D matrix fields to resize
matrixFields = {'c0', 'rho0', 'alpha'};

% Default values for each field
defaultValues = struct('c0', soundSpeedContactMedium, ...
                      'rho0', densityContactMedium, ...
                      'alpha', alphaContactMedium);

% Process each 3D matrix field
for i = 1:length(matrixFields)
    fieldName = matrixFields{i};
    
    if isfield(MaterialMatrix, fieldName)
        currentMatrix = MaterialMatrix.(fieldName);
        [currNx, currNy, currNz] = size(currentMatrix);
        
        % Handle x dimension - beginning first
        if nxPadBegin > 0
            % Pad at the beginning
            padMatrix = defaultValues.(fieldName) * ones(nxPadBegin, currNy, currNz, class(currentMatrix));
            currentMatrix = cat(1, padMatrix, currentMatrix);
            currNx = currNx + nxPadBegin;
        elseif nxPadBegin < 0
            % Trim from the beginning
            currentMatrix = currentMatrix((1 - nxPadBegin):end, :, :);
            currNx = currNx + nxPadBegin;
        end
        
        % Handle x dimension - end
        if nxPadEnd > 0
            % Pad at the end
            padMatrix = defaultValues.(fieldName) * ones(nxPadEnd, currNy, currNz, class(currentMatrix));
            currentMatrix = cat(1, currentMatrix, padMatrix);
            currNx = currNx + nxPadEnd;
        elseif nxPadEnd < 0
            % Trim from the end
            currentMatrix = currentMatrix(1:(currNx + nxPadEnd), :, :);
            currNx = currNx + nxPadEnd;
        end
        
        % Handle y dimension - beginning first
        if nyPadBegin > 0
            % Pad at the beginning
            padMatrix = defaultValues.(fieldName) * ones(currNx, nyPadBegin, currNz, class(currentMatrix));
            currentMatrix = cat(2, padMatrix, currentMatrix);
            currNy = currNy + nyPadBegin;
        elseif nyPadBegin < 0
            % Trim from the beginning
            currentMatrix = currentMatrix(:, (1 - nyPadBegin):end, :);
            currNy = currNy + nyPadBegin;
        end
        
        % Handle y dimension - end
        if nyPadEnd > 0
            % Pad at the end
            padMatrix = defaultValues.(fieldName) * ones(currNx, nyPadEnd, currNz, class(currentMatrix));
            currentMatrix = cat(2, currentMatrix, padMatrix);
            currNy = currNy + nyPadEnd;
        elseif nyPadEnd < 0
            % Trim from the end
            currentMatrix = currentMatrix(:, 1:(currNy + nyPadEnd), :);
            currNy = currNy + nyPadEnd;
        end
        
        % Handle z dimension - beginning first
        if nzPadBegin > 0
            % Pad at the beginning
            padMatrix = defaultValues.(fieldName) * ones(currNx, currNy, nzPadBegin, class(currentMatrix));
            currentMatrix = cat(3, padMatrix, currentMatrix);
            currNz = currNz + nzPadBegin;
        elseif nzPadBegin < 0
            % Trim from the beginning
            currentMatrix = currentMatrix(:, :, (1 - nzPadBegin):end);
            currNz = currNz + nzPadBegin;
        end
        
        % Handle z dimension - end
        if nzPadEnd > 0
            % Pad at the end
            padMatrix = defaultValues.(fieldName) * ones(currNx, currNy, nzPadEnd, class(currentMatrix));
            currentMatrix = cat(3, currentMatrix, padMatrix);
            currNz = currNz + nzPadEnd;
        elseif nzPadEnd < 0
            % Trim from the end
            currentMatrix = currentMatrix(:, :, 1:(currNz + nzPadEnd));
            currNz = currNz + nzPadEnd;
        end
        
        MaterialMatrixResized.(fieldName) = currentMatrix;
    end
end

% Verify final dimensions
if ~isequal(size(MaterialMatrixResized.c0), [nxNew, nyNew, nzNew])
    error('Dimension mismatch after resizing. Expected [%d, %d, %d], got [%d, %d, %d]', ...
        nxNew, nyNew, nzNew, size(MaterialMatrixResized.c0, 1), ...
        size(MaterialMatrixResized.c0, 2), size(MaterialMatrixResized.c0, 3));
end

end

