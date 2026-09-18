function valueString = num2str_significant_figures(value, reference)
%NUM2STR_SIGNIFICANT_FIGURES Convert value to string using reference precision.
%   valueString = num2str_significant_figures(value, reference) returns
%   VALUE formatted with the same number of decimal places as REFERENCE.
%
%   Example:
%       num2str_significant_figures(1.54234, 0.025) returns '1.542'

validateattributes(value, {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, 'value');
validateattributes(reference, {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, 'reference');

% Represent the reference with enough fixed-point precision, then trim
% trailing zeros to recover the effective decimal places.
referenceString = sprintf('%.15f', abs(reference));
referenceString = regexprep(referenceString, '0+$', '');

decimalPointPos = strfind(referenceString, '.');
if isempty(decimalPointPos)
    nDecimalPlaces = 0;
else
    nDecimalPlaces = numel(referenceString) - decimalPointPos;
end

formatString = ['%0.' num2str(nDecimalPlaces) 'f'];
valueString = sprintf(formatString, value);

end
