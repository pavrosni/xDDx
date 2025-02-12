% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [InputStructure] = reshape_transducer_dim(InputStructure)

if isvector(InputStructure.xGrid)
    InputStructure.xGrid = reshape(InputStructure.xGrid,[1 size(InputStructure.xGrid)]);
end

if isvector(InputStructure.yGrid)
    InputStructure.yGrid = reshape(InputStructure.yGrid,[1 size(InputStructure.yGrid)]);
end

if isfield(InputStructure, 'zGrid')
    if isvector(InputStructure.zGrid)
        InputStructure.zGrid = reshape(InputStructure.zGrid,[1 size(InputStructure.zGrid)]);
    end
end

if isfield(InputStructure, 'velocity')
    if ismatrix(InputStructure.velocity)
        InputStructure.velocity = reshape(InputStructure.velocity,[1 size(InputStructure.velocity)]);
    end
end

if isfield(InputStructure, 'pressureWaveforms')
    if ismatrix(InputStructure.pressureWaveforms)
        InputStructure.pressureWaveforms = reshape(InputStructure.pressureWaveforms,[1 size(InputStructure.pressureWaveforms)]);
    end
end

if isfield(InputStructure, 'complexVelocityAmplitude')
    if isvector(InputStructure.complexVelocityAmplitude)
        InputStructure.complexVelocityAmplitude = reshape(InputStructure.complexVelocityAmplitude,[1 size(InputStructure.complexVelocityAmplitude)]);
    end
end

if isfield(InputStructure, 'complexPressureAmplitude')
    if isvector(InputStructure.complexPressureAmplitude)
        InputStructure.complexPressureAmplitude = reshape(InputStructure.complexPressureAmplitude,[1 size(InputStructure.complexPressureAmplitude)]);
    end
end


end

