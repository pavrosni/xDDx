function cudaVersion = normalize_cuda_version(cudaVersion)
%NORMALIZE_CUDA_VERSION Validate CUDA version selection.

if nargin < 1 || isempty(cudaVersion)
    cudaVersion = 'auto';
end

if isnumeric(cudaVersion) && isscalar(cudaVersion)
    cudaVersion = num2str(cudaVersion);
elseif isstring(cudaVersion) && isscalar(cudaVersion)
    cudaVersion = char(cudaVersion);
end

if ~ischar(cudaVersion)
    error('cudaVersion must be ''auto'', ''cuda11'', ''cuda12'', 11, or 12.');
end

cudaVersion = lower(strtrim(cudaVersion));
cudaVersion = regexprep(cudaVersion, '^(cuda|v)?[\s_-]*(11|12)(\.\d+)?$', 'cuda$2');

validCudaVersions = {'auto', 'cuda11', 'cuda12'};
if ~any(strcmp(cudaVersion, validCudaVersions))
    error('cudaVersion must be ''auto'', ''cuda11'', ''cuda12'', 11, or 12.');
end

if ismac && ~strcmp(cudaVersion, 'auto')
    error('cudaVersion = ''%s'' is not supported on macOS. Use CPU mode instead.', cudaVersion);
end

end
