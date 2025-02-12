% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function A = read_matrix_bin(fname)

fid = fopen(fname);
if fid == -1
    error('File is not opened');
end

Ni = fread(fid,1,'int32');
Nj = fread(fid,1,'int32');
Nk = fread(fid,1,'int32');

A = fread(fid,Ni*Nj*Nk,'double');
A = reshape(A,Ni,Nj,Nk);

fclose(fid);

end
