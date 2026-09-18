function options = xddx_backend_options_from_environment()
%XDDX_BACKEND_OPTIONS_FROM_ENVIRONMENT Configure the backend matrix run.

testFolder = fileparts(mfilename('fullpath'));
options = struct();
options.CpuArchitectures = split_selection( ...
    getenv('XDDX_TEST_CPU_ARCHITECTURES'), ...
    {'sse2', 'avx', 'avx2', 'avx512', 'arm64'});
options.CudaVersions = split_selection( ...
    getenv('XDDX_TEST_CUDA_VERSIONS'), {'cuda11', 'cuda12'});
options.ExecutionModes = split_selection( ...
    getenv('XDDX_TEST_EXECUTION_MODES'), {'native', 'docker'});
options.ReportDirectory = strtrim(getenv('XDDX_TEST_REPORT_DIRECTORY'));
if isempty(options.ReportDirectory)
    options.ReportDirectory = fullfile(testFolder, 'reports');
end
end

function values = split_selection(rawValue, defaultValue)
rawValue = strtrim(rawValue);
if isempty(rawValue)
    values = defaultValue;
    return;
end

values = strtrim(strsplit(rawValue, ','));
values = values(~cellfun(@isempty, values));
values = cellfun(@lower, values, 'UniformOutput', false);
end
