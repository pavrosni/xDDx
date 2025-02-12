% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [k, kStep] = make_k_array(numberOfPoints, step)

kStep = 2*pi/(numberOfPoints*step);
if mod(numberOfPoints,2) == 0
    k = ((-numberOfPoints/2):(numberOfPoints/2-1))*kStep;
else
    k = ((-(numberOfPoints-1)/2):((numberOfPoints-1)/2))*kStep;
end

end
