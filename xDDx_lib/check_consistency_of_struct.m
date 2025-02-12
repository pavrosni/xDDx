% Copyright (c) 2025, the code is written by Pavel Rosnitskiy
function isConsistent = check_consistency_of_struct(InputStructure, isTransient)

if isfield(InputStructure,'input')
    if isempty(InputStructure.input)
        InputStructure = rmfield(InputStructure,'input');
    else
        if isTransient
            sizeOfInput = size(InputStructure.input);
            sizeOfInput(end) = [];
            InputStructure.input = zeros([sizeOfInput 1]);
        end
    end
end

if isfield(InputStructure,'xGrid')
    if isvector(InputStructure.xGrid)
       InputStructure.xGrid = reshape(InputStructure.xGrid, [length(InputStructure.xGrid) 1]);
    end
end

if isfield(InputStructure,'yGrid')
    if isvector(InputStructure.yGrid)
       InputStructure.yGrid = reshape(InputStructure.yGrid, [length(InputStructure.yGrid) 1]);
    end
end

if isfield(InputStructure,'zGrid')
    if isvector(InputStructure.zGrid)
       InputStructure.zGrid = reshape(InputStructure.zGrid, [length(InputStructure.zGrid) 1]);
    end
end

if isfield(InputStructure,'input')
    if isvector(InputStructure.input)
       InputStructure.input = reshape(InputStructure.input, [length(InputStructure.input) 1]);
    end
end

if isfield(InputStructure,'dx')
    InputStructure = rmfield(InputStructure,'dx');
end

if isfield(InputStructure,'dy')
    InputStructure = rmfield(InputStructure,'dy');
end

InputStructureCells = struct2cell(InputStructure);
isConsistent = true;

allSizes = [];
for iCells = 1:length(InputStructureCells)
    
    if ~isa(InputStructureCells{iCells},'double')
        isConsistent = false;
    end
    
    try
        sizeInputFull = size(InputStructureCells{iCells});
        allSizes = [allSizes; sizeInputFull];
    catch
        isConsistent = false;
    end
end

if isConsistent
    isConsistent = (nnz(diff(allSizes,1))==0);
end
end

