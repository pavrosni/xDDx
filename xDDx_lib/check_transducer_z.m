% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [Transducer] = check_transducer_z(Transducer)

autoSetZ = true;

if isfield(Transducer, 'radiusOfCurvature')
    if ~isempty(Transducer.radiusOfCurvature)
        radiusOfCurvature = Transducer.radiusOfCurvature;
    end
end

isSphericalSource = exist('radiusOfCurvature', 'var') ~= 0;

if isfield(Transducer, 'zGrid') && (~isempty(Transducer.zGrid))
    autoSetZ = false;
end

if autoSetZ
    if isSphericalSource
        Transducer.zGrid = radiusOfCurvature - sqrt(radiusOfCurvature^2 - Transducer.xGrid.^2 - Transducer.yGrid.^2);
    else
        Transducer.zGrid = zeros(size(Transducer.xGrid));
    end
end

end

