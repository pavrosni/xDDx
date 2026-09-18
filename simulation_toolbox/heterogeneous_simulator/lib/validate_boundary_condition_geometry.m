function isValid = validate_boundary_condition_geometry(MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, ...
    aperture, radiusOfCurvature, radialReserveX, radialReserveY, xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd, useGUI)
%VALIDATE_BOUNDARY_CONDITION_GEOMETRY Validate boundary-condition geometry.
%   Returns true when checks pass. If useGUI is true, shows the setup window
%   and returns false on error. If useGUI is false, throws an error.

    isValid = true;

    [hasError, errMessage] = get_first_error( ...
        MaterialMatrix, izTarget, izBoundaryCondition, radiusOfCurvature, zFieldEnd);

    if ~hasError
        return;
    end

    if useGUI
        show_setup_error_window(errMessage, ...
            MaterialMatrix, ixTarget, iyTarget, izTarget, izBoundaryCondition, ...
            aperture, radiusOfCurvature, radialReserveX, radialReserveY, xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd);
        isValid = false;
    else
        error(errMessage);
    end
end

function [hasError, errMessage] = get_first_error(MaterialMatrix, izTarget, izBoundaryCondition, radiusOfCurvature, zFieldEnd)
    hasError = false;
    errMessage = '';

    isSphericalSource = ~isempty(radiusOfCurvature);

    diffC0 = abs(MaterialMatrix.c0(:, :, izBoundaryCondition) - MaterialMatrix.c0(1, 1, izBoundaryCondition));
    diffRho0 = abs(MaterialMatrix.rho0(:, :, izBoundaryCondition) - MaterialMatrix.rho0(1, 1, izBoundaryCondition));
    if any(diffC0(:) > eps('single')) || any(diffRho0(:) > eps('single'))
        hasError = true;
        errMessage = 'The boundary condition is defined in a non-uniform medium!';
        return;
    end

    if isSphericalSource && (izTarget - izBoundaryCondition) * MaterialMatrix.dz > radiusOfCurvature
        hasError = true;
        errMessage = 'The boundary condition is more distant from the target than the radius of curvature!';
        return;
    end

    if izTarget < izBoundaryCondition
        hasError = true;
        errMessage = 'The target is behind the boundary condition!';
        return;
    end

    if abs(izBoundaryCondition - izTarget) < eps('single')
        hasError = true;
        errMessage = 'The boundary condition is at the target!';
        return;
    end

    if isSphericalSource && zFieldEnd < radiusOfCurvature
        hasError = true;
        errMessage = 'The field end is behind the radius of curvature!';
    end
end
