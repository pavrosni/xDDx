function cases = xddx_backend_matrix(options)
%XDDX_BACKEND_MATRIX Build CPU/GPU variant and launch-mode combinations.

validate_selections(options.CpuArchitectures, ...
    {'sse2', 'avx', 'avx2', 'avx512', 'arm64'}, 'CPU architecture');
validate_selections(options.CudaVersions, ...
    {'cuda11', 'cuda12'}, 'CUDA version');
validate_selections(options.ExecutionModes, ...
    {'native', 'docker'}, 'execution mode');

caseCount = numel(options.ExecutionModes) * ...
    (numel(options.CpuArchitectures) + numel(options.CudaVersions));
emptyCase = struct('Name', '', 'Device', '', 'Variant', '', ...
    'ExecutionMode', '');
cases = repmat(emptyCase, 1, caseCount);
caseIndex = 0;
for modeIndex = 1:numel(options.ExecutionModes)
    mode = options.ExecutionModes{modeIndex};
    for variantIndex = 1:numel(options.CpuArchitectures)
        caseIndex = caseIndex + 1;
        cases(caseIndex) = make_case( ...
            'cpu', options.CpuArchitectures{variantIndex}, mode);
    end
    for variantIndex = 1:numel(options.CudaVersions)
        caseIndex = caseIndex + 1;
        cases(caseIndex) = make_case( ...
            'cuda', options.CudaVersions{variantIndex}, mode);
    end
end
end

function value = make_case(device, variant, mode)
value = struct('Name', [mode '_' device '_' variant], ...
    'Device', device, 'Variant', variant, 'ExecutionMode', mode);
end

function validate_selections(values, allowedValues, label)
if ~iscell(values) || any(~ismember(values, allowedValues))
    error('xDDx:Tests:InvalidBackendSelection', ...
        'Invalid %s selection. Allowed values: %s.', ...
        label, strjoin(allowedValues, ', '));
end
end
