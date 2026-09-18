function [prepareSimulationOnly, presimulatedOutputFile, preparedSimulationDataPath] = ...
    plan_presimulation_workflow(kWaveCalculationFlag, prepareSimulationOnly, showPresimulatedData, preparedSimulationDataPath, ~)
%PLAN_PRESIMULATION_WORKFLOW Validate and normalize the CPU presimulation workflow options.

if prepareSimulationOnly && ~strcmpi(kWaveCalculationFlag, 'cpu')
    error('prepareSimulationOnly is supported only when kWaveCalculationFlag = ''cpu''.');
end

presimulatedOutputFile = '';
if ~isequal(showPresimulatedData, false)
    if ~strcmpi(kWaveCalculationFlag, 'cpu')
        error('showPresimulatedData is supported only when kWaveCalculationFlag = ''cpu''.');
    end
    presimulatedOutputFile = char(showPresimulatedData);
end

if prepareSimulationOnly && ~isempty(presimulatedOutputFile)
    error('prepareSimulationOnly and showPresimulatedData cannot be used at the same time.');
end

end
