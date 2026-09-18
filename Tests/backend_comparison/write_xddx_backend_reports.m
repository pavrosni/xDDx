function reportPaths = write_xddx_backend_reports( ...
    results, hostInfo, options, runInfo)
%WRITE_XDDX_BACKEND_REPORTS Write incremental CSV and MAT reports.

if exist(options.ReportDirectory, 'dir') ~= 7
    mkdir(options.ReportDirectory);
end

baseName = ['xddx_backend_comparison_' runInfo.RunId];
reportPaths = struct( ...
    'Csv', fullfile(options.ReportDirectory, [baseName '.csv']), ...
    'Mat', fullfile(options.ReportDirectory, [baseName '.mat']));

reportTable = struct2table(results);
rowCount = height(reportTable);
reportTable.RunId = repmat(string(runInfo.RunId), rowCount, 1);
reportTable.StartedUtc = repmat(string(runInfo.StartedUtc), rowCount, 1);
reportTable.ReportSchemaVersion = repmat( ...
    runInfo.ReportSchemaVersion, rowCount, 1);
reportTable.CaseDefinition = repmat( ...
    string(runInfo.CaseDefinition), rowCount, 1);
reportTable.MachineName = repmat(string(hostInfo.MachineName), rowCount, 1);
reportTable.Platform = repmat(string(hostInfo.Platform), rowCount, 1);
reportTable.Architecture = repmat(string(hostInfo.Architecture), rowCount, 1);
reportTable.OperatingSystem = repmat( ...
    string(hostInfo.OperatingSystem), rowCount, 1);
reportTable.IsWsl = repmat(hostInfo.IsWsl, rowCount, 1);
reportTable.CpuModel = repmat(string(hostInfo.CpuModel), rowCount, 1);
reportTable.NvidiaSummary = repmat( ...
    string(hostInfo.NvidiaSummary), rowCount, 1);
reportTable.DockerVersion = repmat( ...
    string(hostInfo.DockerVersion), rowCount, 1);
reportTable.DockerAvailable = repmat( ...
    hostInfo.DockerAvailable, rowCount, 1);
reportTable.MaximumCudaVersion = repmat( ...
    hostInfo.MaximumCudaVersion, rowCount, 1);
reportTable.MatlabVersion = repmat( ...
    string(hostInfo.MatlabVersion), rowCount, 1);
reportTable.MatlabRelease = repmat( ...
    string(hostInfo.MatlabRelease), rowCount, 1);
reportTable.GitCommit = repmat(string(hostInfo.GitCommit), rowCount, 1);

writetable(reportTable, reportPaths.Csv);
save(reportPaths.Mat, 'results', 'hostInfo', 'options', ...
    'runInfo', 'reportPaths', '-v7.3');
end
