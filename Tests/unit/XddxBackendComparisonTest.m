classdef XddxBackendComparisonTest < matlab.unittest.TestCase
    %XDDXBACKENDCOMPARISONTEST Backend matrix and report utilities.

    methods (TestClassSetup)
        function addSupportPaths(testCase)
            testRoot = fileparts(fileparts(mfilename('fullpath')));
            projectRoot = fileparts(testRoot);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testRoot, 'backend_comparison')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'xDDx_lib')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'simulation_toolbox', ...
                ['heterogeneous_' 'simulator'], 'lib')));
        end
    end

    methods (Test)
        function testCompleteBackendMatrixHasFourteenCases(testCase)
            options = complete_options();

            actual = xddx_backend_matrix(options);

            testCase.verifyNumElements(actual, 14);
            testCase.verifyEqual({actual.Name}, { ...
                'native_cpu_sse2', 'native_cpu_avx', ...
                'native_cpu_avx2', 'native_cpu_avx512', ...
                'native_cpu_arm64', 'native_cuda_cuda11', ...
                'native_cuda_cuda12', 'docker_cpu_sse2', ...
                'docker_cpu_avx', 'docker_cpu_avx2', ...
                'docker_cpu_avx512', 'docker_cpu_arm64', ...
                'docker_cuda_cuda11', 'docker_cuda_cuda12'});
        end

        function testOutputAndDifferenceNorms(testCase)
            output = [3, 4];
            reference = [0, 4];

            actual = xddx_backend_norms(output, reference);

            testCase.verifyEqual(actual.OutputL2, 5, 'AbsTol', eps);
            testCase.verifyEqual(actual.OutputLInf, 4, 'AbsTol', eps);
            testCase.verifyEqual(actual.DifferenceL2, 3, 'AbsTol', eps);
            testCase.verifyEqual(actual.DifferenceLInf, 3, 'AbsTol', eps);
            testCase.verifyEqual(actual.RelativeL2, 0.75, 'AbsTol', eps);
            testCase.verifyEqual(actual.RelativeLInf, 0.75, 'AbsTol', eps);
        end

        function testSimulatorOverridesChangeOnlyRequestedFields(testCase)
            defaults = struct('kWaveCalculationFlag', 'auto', ...
                'xDDxCalculationFlag', 'auto', 'CFL', 0.3, ...
                'useGUI', true, 'shouldPlot', true);
            overrides = struct('kWaveCalculationFlag', 'cpu', ...
                'xDDxCalculationFlag', 'cpu', ...
                'useGUI', false, 'shouldPlot', false);

            actual = apply_xddx_simulator_input_overrides( ...
                defaults, overrides);

            testCase.verifyEqual(actual.CFL, 0.3, 'AbsTol', eps);
            testCase.verifyEqual(actual.kWaveCalculationFlag, 'cpu');
            testCase.verifyEqual(actual.xDDxCalculationFlag, 'cpu');
            testCase.verifyFalse(actual.useGUI);
            testCase.verifyFalse(actual.shouldPlot);
        end

        function testUnknownSimulatorOverrideIsRejected(testCase)
            defaults = struct('CFL', 0.3);
            overrides = struct('CFL', 0.1);

            testCase.verifyError(@() ...
                apply_xddx_simulator_input_overrides(defaults, overrides), ...
                'xDDx:Simulator:InvalidOverrideField');
        end

        function testUnavailableDockerCaseReturnsSkip(testCase)
            backendCase = struct('Name', 'docker_cpu_sse2', ...
                'Device', 'cpu', 'Variant', 'sse2', ...
                'ExecutionMode', 'docker');
            hostInfo = struct('DockerAvailable', false, ...
                'DockerMessage', 'Docker is not installed.');

            [available, reason] = xddx_backend_case_availability( ...
                backendCase, hostInfo, '');

            testCase.verifyFalse(available);
            testCase.verifySubstring(reason, 'Docker is unavailable');
        end

        function testDiscoveredDockerCliDirectoryIsAddedToPath(testCase)
            temporaryFolder = string(tempname);
            mkdir(temporaryFolder);
            testCase.addTeardown(@() rmdir(temporaryFolder, 's'));
            dockerExecutable = fullfile(temporaryFolder, 'docker');
            fileId = fopen(dockerExecutable, 'w');
            testCase.assertNotEqual(fileId, -1);
            fclose(fileId);
            originalPath = getenv('PATH');
            testCase.addTeardown(@() setenv('PATH', originalPath));

            fix = ensure_xddx_docker_cli_in_path({dockerExecutable});

            testCase.verifyTrue(fix.didChangePath);
            testCase.verifyEqual(fix.addedDirs, {char(temporaryFolder)});
            pathDirectories = regexp(getenv('PATH'), pathsep, 'split');
            testCase.verifyEqual(pathDirectories{1}, char(temporaryFolder));
        end
    end
end

function options = complete_options()
options = struct();
options.CpuArchitectures = ...
    {'sse2', 'avx', 'avx2', 'avx512', 'arm64'};
options.CudaVersions = {'cuda11', 'cuda12'};
options.ExecutionModes = {'native', 'docker'};
end
