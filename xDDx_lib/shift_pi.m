% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [out] = shift_pi(in)

notWithinPi = (abs(in) > pi);

inShift = in(notWithinPi) + pi;

isPositive = (inShift > 0);

inShift = mod(inShift, 2*pi);
inShift(isPositive & (inShift==0)) = 2*pi;

out = in;
out(notWithinPi) = inShift - pi;

end

