function metrics = xddx_backend_norms(output, reference)
%XDDX_BACKEND_NORMS Calculate output and reference-difference norms.

if ~isnumeric(output) || isempty(output) || any(~isfinite(output(:)))
    error('xDDx:Tests:InvalidBackendOutput', ...
        'The simulator output must be a nonempty finite numeric array.');
end
if nargin < 2
    reference = [];
end

metrics = struct();
metrics.OutputL2 = double(norm(output(:), 2));
metrics.OutputLInf = double(norm(output(:), Inf));
metrics.DifferenceL2 = NaN;
metrics.DifferenceLInf = NaN;
metrics.RelativeL2 = NaN;
metrics.RelativeLInf = NaN;

if isempty(reference)
    return;
end
if ~isequal(size(output), size(reference))
    error('xDDx:Tests:BackendOutputSizeMismatch', ...
        'Backend output size differs from the reference output size.');
end

difference = output(:) - reference(:);
referenceL2 = double(norm(reference(:), 2));
referenceLInf = double(norm(reference(:), Inf));
metrics.DifferenceL2 = double(norm(difference, 2));
metrics.DifferenceLInf = double(norm(difference, Inf));
metrics.RelativeL2 = metrics.DifferenceL2 / max(referenceL2, eps);
metrics.RelativeLInf = metrics.DifferenceLInf / max(referenceLInf, eps);
end
