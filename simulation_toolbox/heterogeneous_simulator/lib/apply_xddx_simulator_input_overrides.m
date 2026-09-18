function simulatorInputs = apply_xddx_simulator_input_overrides( ...
    simulatorInputs, overrides)
%APPLY_XDDX_SIMULATOR_INPUT_OVERRIDES Apply validated script overrides.
%   This keeps xDDx_simulator.m as the source of truth for its default case
%   while allowing automated runs to select a backend and disable UI output.

if isempty(overrides)
    return;
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('xDDx:Simulator:InvalidOverrides', ...
        'simulatorInputOverrides must be a scalar structure.');
end

allowedFields = {'kWaveCalculationFlag', 'xDDxCalculationFlag', ...
    'cpuArchitecture', 'cudaVersion', 'useGUI', 'shouldPlot'};
overrideFields = fieldnames(overrides);
unknownFields = setdiff(overrideFields, allowedFields);
if ~isempty(unknownFields)
    error('xDDx:Simulator:InvalidOverrideField', ...
        'Unsupported simulator input override(s): %s.', ...
        strjoin(unknownFields, ', '));
end

for fieldIndex = 1:numel(overrideFields)
    fieldName = overrideFields{fieldIndex};
    simulatorInputs.(fieldName) = overrides.(fieldName);
end
end
