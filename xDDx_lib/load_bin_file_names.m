% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [BinFileNames] = load_bin_file_names

BinFileNames = [];
BinFileNames.xSource    = 'x_source.bin';
BinFileNames.ySource    = 'y_source.bin';
BinFileNames.zSource    = 'z_source.bin';
BinFileNames.xField    = 'x_field.bin';
BinFileNames.yField    = 'y_field.bin';
BinFileNames.zField    = 'z_field.bin';
BinFileNames.reInput  = 're_input.bin';
BinFileNames.imInput  = 'im_input.bin';
BinFileNames.reOutput = 're_output.bin';
BinFileNames.imOutput = 'im_output.bin';
BinFileNames.vectorParam = 'vec_input_param.bin';

end

