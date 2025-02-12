% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [out] = shift_2pi(in)

isPositive = (in > 0);

out = in;

out = mod(out, 2*pi);
out(isPositive & (out==0)) = 2*pi;


end

