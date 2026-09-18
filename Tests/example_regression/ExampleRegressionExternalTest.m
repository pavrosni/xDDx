classdef ExampleRegressionExternalTest < matlab.unittest.TestCase
    %EXAMPLEREGRESSIONEXTERNALTEST Optional old-vs-reorganized comparison.

    methods (TestClassSetup)
        function addSupportPath(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fileparts(mfilename('fullpath'))));
        end
    end

    methods (Test, TestTags={'External', 'ExampleRegression'})
        function testSelectedExamplesMatchLegacyProject(testCase)
            selection = getenv('XDDX_EXAMPLE_REGRESSION_SCRIPTS');
            testCase.assumeNotEmpty(selection, ...
                ['Set XDDX_EXAMPLE_REGRESSION_SCRIPTS to a comma-separated ' ...
                'list of mapped example names before running external tests.']);
            referenceProject = string(getenv( ...
                'XDDX_EXAMPLE_REGRESSION_REFERENCE'));
            referenceProject(strlength(referenceProject) == 0) = string('xddx');

            [summary, comparison, ~, results] = run_example_regression_comparison( ...
                'Scripts', strtrim(strsplit(selection, ',')), ...
                'SimulationDevices', 'cpu', ...
                'ReferenceProject', referenceProject, ...
                'SaveReports', true);

            testCase.verifyTrue(all(summary.Success));
            testCase.verifyEqual({results.referenceProject}, ...
                repmat({char(referenceProject)}, 1, numel(results)));
            comparable = comparison.SameSize;
            testCase.verifyTrue(any(comparable));
            testCase.verifyLessThanOrEqual( ...
                max(comparison.RelativeL2(comparable)), 1e-9);
        end
    end
end
