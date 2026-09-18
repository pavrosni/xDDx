classdef TransientBackendExternalTest < matlab.unittest.TestCase
    %TRANSIENTBACKENDEXTERNALTEST Optional native/Docker backend comparison.

    methods (TestClassSetup)
        function addTestSupportPath(testCase)
            testRoot = fileparts(fileparts(mfilename('fullpath')));
            projectRoot = fileparts(testRoot);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'xDDx_lib')));
        end
    end

    methods (Test, TestTags={'External', 'BackendComparison'})
        function testRequestedTransientBackendsAgree(testCase)
            options = transient_backend_options_from_environment();
            results = run_transient_backend_comparison(options);
            successful = results([results.Success]);

            testCase.verifyNotEmpty(successful, ...
                'No requested transient backend completed successfully.');
            testCase.verifyTrue(all([results.Success]), ...
                'At least one requested transient backend failed.');
            testCase.verifyLessThanOrEqual(max([successful.RelativeL2]), ...
                options.RelativeTolerance);
        end
    end
end
