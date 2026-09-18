function [nx, ny, nz] = adjust_fft_size(nx, ny, nz, pml_x_size, pml_y_size, pml_z_size)
%ADJUST_FFT_SIZE Make 3-D FFT grid sizes even and FFT-friendly.
%   [nx, ny, nz] = ADJUST_FFT_SIZE(nx, ny, nz) returns the sizes
%   greater than or equal to the requested sizes such that nx, ny, and nz are even and factors only into {2, 3, 5}.
nx = nx + 2*pml_x_size;
ny = ny + 2*pml_y_size;
nz = nz + 2*pml_z_size;

% The FFT grid including the PML must be both even and factorable by 2, 3, and 5.
if mod(nx, 2) ~= 0 || ~isFastLength(nx)
nx = nextFastEven(nx);
end
if mod(ny, 2) ~= 0 || ~isFastLength(ny)
    ny = nextFastEven(ny);
end
if mod(nz, 2) ~= 0 || ~isFastLength(nz)
    nz = nextFastEven(nz);
end

nx = nx - 2*pml_x_size;
ny = ny - 2*pml_y_size;
nz = nz - 2*pml_z_size;

end

function n = nextFastEven(n)
%NEXTFASTEVEVN Return the next even length with factors in {2,3,5}.
    validateattributes(n, {'numeric'}, {'scalar', 'integer', '>=', 1});

    % Start at the next even length not smaller than n
    n = double(n);
    if mod(n, 2) ~= 0
        n = n + 1;
    end

    while ~isFastLength(n)
        n = n + 2; %#ok<AGROW> % keep parity even
    end
end

function tf = isFastLength(n)
%ISFASTLENGTH True if n has no prime factors other than 2, 3, or 5.
    if n == 1
        tf = true;
        return;
    end
    for p = [2, 3, 5]
        while mod(n, p) == 0
            n = n / p;
        end
    end
    tf = (n == 1);
end


