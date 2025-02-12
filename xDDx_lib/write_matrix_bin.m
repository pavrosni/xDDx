% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function write_matrix_bin(fname, A)

fid = fopen(fname,'w');
if fid == -1
    error('File is not opened');
end

fwrite(fid,size(A,1),'int32');
fwrite(fid,size(A,2),'int32');
fwrite(fid,size(A,3),'int32');
fwrite(fid,A(:),'double');
fclose(fid);

end
