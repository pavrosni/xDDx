function results = run_transient_backend_comparison(options)
%RUN_TRANSIENT_BACKEND_COMPARISON Compare tiny transient backend outputs.

cases = build_cases(options);
if isempty(cases)
    error('xDDx:Tests:NoBackends', 'No backend cases were selected.');
end

oldDockerValue = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
cleanup = onCleanup(@() setenv('XDDX_USE_DOCKER_ON_WINDOWS', oldDockerValue));

source = transient_source();
field = field_points();
medium = struct('soundSpeed', 1500, 'density', 1000);
results = repmat(empty_result(), 1, numel(cases));

reference = [];
for caseIndex = 1:numel(cases)
    current = cases(caseIndex);
    setenv('XDDX_USE_DOCKER_ON_WINDOWS', char(string(strcmp(current.ExecutionMode, 'docker'))));
    service = struct('threadsPerBlockGPU', 128);
    if strcmp(current.Device, 'cpu')
        service.cpuArchitecture = current.Variant;
    else
        service.cudaVersion = current.Variant;
    end

    results(caseIndex).Name = current.Name;
    results(caseIndex).Device = current.Device;
    results(caseIndex).Variant = current.Variant;
    results(caseIndex).ExecutionMode = current.ExecutionMode;

    try
        timer = tic;
        output = rayleigh_simulator(1, 0.5e6, 4, current.Device, true, ...
            source, field, medium, service);
        results(caseIndex).ElapsedSeconds = toc(timer);
        results(caseIndex).Success = true;
        results(caseIndex).Output = output;
        if isempty(reference)
            reference = output;
        end
        results(caseIndex).RelativeL2 = norm(output(:) - reference(:))/ ...
            max(norm(reference(:)), eps);
    catch exception
        results(caseIndex).ErrorIdentifier = exception.identifier;
        results(caseIndex).ErrorMessage = exception.message;
    end
end
end

function cases = build_cases(options)
cases = struct('Name', {}, 'Device', {}, 'Variant', {}, 'ExecutionMode', {});
for modeIndex = 1:numel(options.ExecutionModes)
    mode = options.ExecutionModes{modeIndex};
    for variantIndex = 1:numel(options.CpuArchitectures)
        cases(end + 1) = make_case('cpu', ...
            options.CpuArchitectures{variantIndex}, mode); %#ok<AGROW>
    end
    for variantIndex = 1:numel(options.CudaVersions)
        cases(end + 1) = make_case('cuda', ...
            options.CudaVersions{variantIndex}, mode); %#ok<AGROW>
    end
end
end

function value = make_case(device, variant, mode)
value = struct('Name', [mode '_' device '_' variant], ...
    'Device', device, 'Variant', variant, 'ExecutionMode', mode);
end

function source = transient_source()
[xGrid, yGrid] = build_flat_grid_centered(3, 1e-3, 3, 1e-3);
baseInput = ones(size(xGrid))/(1500*1000);
source = struct('xGrid', xGrid, 'yGrid', yGrid, ...
    'zGrid', zeros(size(xGrid)), 'dx', 1e-3, 'dy', 1e-3, ...
    'input', cat(3, baseInput, 0.5*baseInput));
end

function field = field_points()
field = struct('xGrid', reshape([0, 1e-3], 1, 1, []), ...
    'yGrid', zeros(1, 1, 2), ...
    'zGrid', reshape([0.02, 0.03], 1, 1, []));
end

function value = empty_result()
value = struct('Name', '', 'Device', '', 'Variant', '', ...
    'ExecutionMode', '', 'Success', false, 'ElapsedSeconds', NaN, ...
    'RelativeL2', NaN, 'ErrorIdentifier', '', 'ErrorMessage', '', ...
    'Output', []);
end
