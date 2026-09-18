classdef XddxBackendExternalTest < matlab.unittest.TestCase
    %XDDXBACKENDEXTERNALTEST Default-case backend comparison on this PC.

    methods (TestClassSetup)
        function addTestPaths(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            projectRoot = fileparts(fileparts(testFolder));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                testFolder));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'xDDx_lib')));
        end
    end

    methods (Test, TestTags={'External', 'BackendComparison', 'Slow'})
        function testDefaultCaseAcrossAvailableBackends(testCase)
            options = xddx_backend_options_from_environment();
            [results, reportPaths] = run_xddx_backend_comparison(options);
            attempted = results([results.Attempted]);
            successful = results([results.Success]);

            testCase.assumeNotEmpty(attempted, ...
                'No selected backend is runnable on this test PC.');
            testCase.verifyTrue(all([attempted.Success]), ...
                'At least one runnable backend failed; inspect the report.');
            testCase.verifyNotEmpty(successful);
            testCase.verifyTrue(all(isfinite([successful.OutputL2])));
            testCase.verifyTrue(all(isfinite([successful.OutputLInf])));
            testCase.verifyEqual(exist(reportPaths.Csv, 'file'), 2);
            testCase.verifyEqual(exist(reportPaths.Mat, 'file'), 2);
        end
    end
end
