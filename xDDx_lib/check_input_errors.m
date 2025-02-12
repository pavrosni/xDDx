% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [inputRayleigh, surfElementArea, isInputSource, isInputField] = check_input_errors(simulationDevice, isTransient, regime, regimesAllVector, errorMessages,SourceParameters,FieldParameters,spatialTolerance, radiusOfCurvature)

isInputSource = false;
isInputField = false;

if ~any(strcmpi(simulationDevice,{'cpu', 'cuda'}))
    error(errorMessages.simulationDevice);
end

if ~ismember(regime,regimesAllVector)
    error(errorMessages.regime);
end

try
    inputSource = SourceParameters.input;
catch
    inputSource = [];
end

try
    inputField = FieldParameters.input;
catch
    inputField = [];
end

if ~check_consistency_of_struct(SourceParameters, isTransient)
  error(errorMessages.checkSourceParametersConsistency);   
end

if ~check_consistency_of_struct(FieldParameters, isTransient)
  error(errorMessages.checkFieldParametersConsistency);   
end

if max(SourceParameters.zGrid(:)) >= min(FieldParameters.zGrid(:))
    error(errorMessages.sourcePosition);
end

switch regime
    case 1 % 1 Backpropagation: P on a plane --> V on a plane
        
        if ~check_plane_fit(FieldParameters.zGrid, spatialTolerance)
            error(errorMessages.fieldGridPlane);
        end
        if ~check_plane_fit(SourceParameters.zGrid, spatialTolerance)
            error(errorMessages.sourceGridPlane);
        end
        if (isempty(inputField) || ~isempty(inputSource))
            error(errorMessages.checkInputStructure);
        end
        
        
    case 2 % 2 Backpropagation: P on a plane --> V on a sphere
        if max(SourceParameters.xGrid(:).^2+SourceParameters.yGrid(:).^2) > radiusOfCurvature^2
           error(errorMessages.sourceCurvature); 
        end
        if ~check_plane_fit(FieldParameters.zGrid, spatialTolerance)
            error(errorMessages.fieldGridPlane);
        end
        if ~check_sphere_fit(SourceParameters.xGrid, SourceParameters.yGrid, SourceParameters.zGrid, spatialTolerance, radiusOfCurvature)
            error(errorMessages.sourceGridSphere);
        end
        if (isempty(inputField) || ~isempty(inputSource))
            error(errorMessages.checkInputStructure);
        end
        
    case 3 % 3 Forwardpropagation: V on a sphere --> P at an arbitrary set of points
        if max(SourceParameters.xGrid(:))^2 + max(SourceParameters.yGrid(:))^2 > radiusOfCurvature^2
           error(errorMessages.sourceCurvature); 
        end
        
        if ~check_sphere_fit(SourceParameters.xGrid, SourceParameters.yGrid, SourceParameters.zGrid, spatialTolerance, radiusOfCurvature)
            error(errorMessages.sourceGridSphere);
        end
        if (~isempty(inputField) || isempty(inputSource))
            error(errorMessages.checkInputStructure);
        end
        
    case 4 %CalcPlaneSourceInitialV% 4 Forwardpropagation: V on a plane --> P at an arbitrary set of points
        
        if ~check_plane_fit(SourceParameters.zGrid, spatialTolerance)
            error(errorMessages.sourceGridPlane);
        end
        if (~isempty(inputField) || isempty(inputSource))
            error(errorMessages.checkInputStructure);
        end
        
    case 5 % 5 Forwardpropagation: P on a plane --> P at an arbitrary set of points
        
        if ~check_plane_fit(SourceParameters.zGrid, spatialTolerance)
            error(errorMessages.sourceGridPlane);
        end
        if (~isempty(inputField) || isempty(inputSource))
            error(errorMessages.checkInputStructure);
        end
        
    case 6 % 6 Backpropagation: V on a sphere --> P on a planes
        if max(FieldParameters.xGrid(:))^2+max(FieldParameters.yGrid(:))^2 > radiusOfCurvature^2
           error(errorMessages.fieldCurvature); 
        end
        
        if ~check_sphere_fit(FieldParameters.xGrid, FieldParameters.yGrid, FieldParameters.zGrid, spatialTolerance, radiusOfCurvature)
            error(errorMessages.fieldGridSphere);
        end
        if ~check_plane_fit(SourceParameters.zGrid, spatialTolerance)
            error(errorMessages.sourceGridPlane);
        end
        if (isempty(inputField) || ~isempty(inputSource))
            error(errorMessages.checkInputStructure);
        end
end


if ~isempty(inputSource)
    inputRayleigh = SourceParameters.input;
    surfElementArea = SourceParameters.dx*SourceParameters.dy;
    isInputSource = true;
elseif ~isempty(inputField)
    inputRayleigh = FieldParameters.input;
    surfElementArea = FieldParameters.dx*FieldParameters.dy;
    isInputField = true;
else
    error(errorMessages.noInput);
end
 
    

end

