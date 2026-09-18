function validate_field_limits(xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, zFieldBegin, zFieldEnd)
%VALIDATE_FIELD_LIMITS Check that Begin <= End for x, y, and z field limits
%   VALIDATE_FIELD_LIMITS(xFieldBegin, xFieldEnd, yFieldBegin, yFieldEnd, ...
%                         zFieldBegin, zFieldEnd)
%   Validates that each Begin value is less than or equal to its 
%   corresponding End value. Throws an error with a descriptive message 
%   if any validation fails.
%
%   Example:
%       validate_field_limits(-15e-3, 15e-3, 0, 0, 40e-3, 80e-3)

    % Check x dimension
    if xFieldBegin > xFieldEnd
        error('validate_field_limits:InvalidX', ...
              'xFieldBegin (%.6e) must be <= xFieldEnd (%.6e)', ...
              xFieldBegin, xFieldEnd);
    end
    
    % Check y dimension
    if yFieldBegin > yFieldEnd
        error('validate_field_limits:InvalidY', ...
              'yFieldBegin (%.6e) must be <= yFieldEnd (%.6e)', ...
              yFieldBegin, yFieldEnd);
    end
    
    % Check z dimension
    if zFieldBegin > zFieldEnd
        error('validate_field_limits:InvalidZ', ...
              'zFieldBegin (%.6e) must be <= zFieldEnd (%.6e)', ...
              zFieldBegin, zFieldEnd);
    end
    
end

