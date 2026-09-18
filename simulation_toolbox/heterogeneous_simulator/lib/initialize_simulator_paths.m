function kWaveBinPath = initialize_simulator_paths(scriptDirectory, kWavePath)
%INITIALIZE_SIMULATOR_PATHS Resolve and add simulator library paths.

kWavePath = char(kWavePath);
if java.io.File(kWavePath).isAbsolute()
    kWaveRoot = char(java.io.File(kWavePath).getCanonicalPath());
else
    kWaveRoot = resolve_script_relative_path(scriptDirectory, kWavePath);
end
kWaveBinPath = fullfile(kWaveRoot, 'binaries');

addpath(genpath(kWaveRoot));

end
