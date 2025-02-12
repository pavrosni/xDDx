% Copyright (c) 2025, the code is written by Pavel Rosnitskiy

function [errorMessages] = load_error_messages

errorMessages = [];
errorMessages.simulationDevice = 'simulationDevice must be either "cuda" or "cpu"!';
errorMessages.sourceGridSphere = 'Please check SourceParameters! It seems that the grid points are not on a sphere, or radiusOfCurvature does not fit the sphere.';
errorMessages.sourceGridPlane = 'Please check SourceParameters! It seems that the grid points are not on a plane.';
errorMessages.fieldGridSphere = 'Please check FieldParameters! It seems that the grid points are not on a sphere, or radiusOfCurvature does not fit the sphere.';
errorMessages.fieldGridPlane = 'Please check FieldParameters! It seems that the grid points are not on a plane.';
errorMessages.regime = 'Wrong simulation regime!';
errorMessages.exe = 'Rayleigh Integral Simulator cannot be executed!';
errorMessages.noInput = 'No inputField!';
errorMessages.checkInputStructure = "Check SourceParameters.input or FieldParameters.input!";
errorMessages.checkSourceParametersConsistency = 'The fields of SourceParameters structure are of inconsistent size or non-double data type!';
errorMessages.checkFieldParametersConsistency = 'The fields of FieldParameters structure are of inconsistent size or non-double data type!';
errorMessages.sourcePosition = 'The z-position of the Source cannot be greater than that of the Field!';
errorMessages.sourceCurvature = 'The on-sphere grid for Source cannot be generated for a given focalLength! Please check focalLength, SourceParameters.xGrid, or SourceParameters.yGrid!';
errorMessages.fieldCurvature = 'The on-sphere grid for Field cannot be generated for a given focalLength! Please check focalLength, FieldParameters.xGrid, or FieldParameters.yGrid!';
errorMessages.noCuda = 'The simulation has not been performed! The most likely issue is that your graphics card is not CUDA-compatible. Try running your simulation in CPU mode (simulationDevice = ''cpu''). If you are sure that your GPU is CUDA-compatible, please check that you have the latest version of your graphics drivers.';
errorMessages.noLinux = 'Linux OS is not supported yet!';
errorMessages.noMac = 'Mac OS is not supported yet!';
errorMessages.noOS = 'Your operating system is not supported yet!';
errorMessages.simInput = 'Check input parameters for rayleigh_simulator!';

end

