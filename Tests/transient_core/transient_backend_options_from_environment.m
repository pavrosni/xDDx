function options = transient_backend_options_from_environment()
%TRANSIENT_BACKEND_OPTIONS_FROM_ENVIRONMENT Configure optional backend tests.

options = struct();
options.CpuArchitectures = split_selection( ...
    getenv('XDDX_TEST_CPU_ARCHITECTURES'), {'avx2'});
options.CudaVersions = split_selection( ...
    getenv('XDDX_TEST_CUDA_VERSIONS'), {});
options.ExecutionModes = split_selection( ...
    getenv('XDDX_TEST_EXECUTION_MODES'), {'native'});
options.RelativeTolerance = 1e-9;
end

function values = split_selection(rawValue, defaultValue)
rawValue = strtrim(rawValue);
if isempty(rawValue)
    values = defaultValue;
    return;
end

values = strtrim(strsplit(rawValue, ','));
values = values(~cellfun(@isempty, values));
end
