function results = run_all_tests(varargin)
%RUN_ALL_TESTS Run the xDDx test suite, excluding external integrations by default.
%
%   results = run_all_tests
%   results = run_all_tests('IncludeExternal', true)
%
% External tests include Docker execution and full legacy example runs. They
% remain discoverable, but the default run excludes them because they require
% services or can take a long time. Heterogeneous simulator source folders are
% not included in test discovery; its external backend test invokes the public
% entry script explicitly. Tests run serially because native binaries use fixed
% temporary filenames.

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'IncludeExternal', false, ...
    @(value) islogical(value) && isscalar(value));
parse(parser, varargin{:});

testRoot = fileparts(mfilename('fullpath'));

if parser.Results.IncludeExternal
    results = runtests(testRoot, 'IncludeSubfolders', true, ...
        'UseParallel', false);
else
    normalTestFolders = {fullfile(testRoot, 'unit'), ...
        fullfile(testRoot, 'integration')};
    results = matlab.unittest.TestResult.empty;
    for folderIndex = 1:numel(normalTestFolders)
        folderResults = runtests(normalTestFolders{folderIndex}, ...
            'IncludeSubfolders', true, 'UseParallel', false);
        results = [results, folderResults]; %#ok<AGROW>
    end
end

disp(table(results));
end
