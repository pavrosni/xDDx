classdef XddxIoAndConfigurationTest < matlab.unittest.TestCase
    %XDDXIOANDCONFIGURATIONTEST File formats and backend configuration.

    properties (SetAccess=private)
        ProjectRoot
    end

    methods (TestClassSetup)
        function addLibraryPath(testCase)
            testCase.ProjectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            libraryDir = fullfile(testCase.ProjectRoot, 'xDDx_lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(libraryDir));
        end
    end

    methods (Test)
        function testBinaryMatrixRoundTrip(testCase)
            temporaryFolder = string(tempname);
            mkdir(temporaryFolder);
            testCase.addTeardown(@() rmdir(temporaryFolder, 's'));
            filePath = fullfile(temporaryFolder, 'matrix.bin');
            expected = reshape(1:24, [2, 3, 4]);

            write_matrix_bin(filePath, expected);
            actual = read_matrix_bin(filePath);

            testCase.verifyEqual(actual, expected, 'AbsTol', eps);
        end

        function testBinaryFileNamesAreUnique(testCase)
            actual = load_bin_file_names();
            fileNames = struct2cell(actual);

            testCase.verifyEqual(numel(unique(fileNames)), numel(fileNames));
            testCase.verifyTrue(all(endsWith(fileNames, '.bin')));
        end

        function testErrorMessageCatalogContainsSimulatorFailures(testCase)
            actual = load_error_messages();

            testCase.verifyTrue(isfield(actual, 'simulationDevice'));
            testCase.verifyTrue(isfield(actual, 'dockerRun'));
            testCase.verifyTrue(isfield(actual, 'noCuda'));
        end

        function testSaveTransducerWritesExpectedStructure(testCase)
            temporaryFolder = string(tempname);
            mkdir(temporaryFolder);
            testCase.addTeardown(@() rmdir(temporaryFolder, 's'));
            oldFolder = cd(temporaryFolder);
            testCase.addTeardown(@() cd(oldFolder));
            xGrid = [0, 1e-3];
            yGrid = [0, 0];
            zGrid = [0, 0];
            velocity = [1, 2];

            save_transducer(false, 1, 1e6, xGrid, yGrid, zGrid, ...
                [], 1e-3, 1e-3, velocity);
            files = dir(fullfile(temporaryFolder, 'transducer_*.mat'));
            loaded = load(fullfile(files(1).folder, files(1).name), 'TransducerSf');

            testCase.verifyEqual(numel(files), 1);
            testCase.verifyEqual(loaded.TransducerSf.complexVelocityAmplitude, ...
                velocity, 'AbsTol', eps);
            testCase.verifyFalse(isfield(loaded.TransducerSf, 'radiusOfCurvature'));
        end

        function testTransducerSpreadsheetTemplateLoads(testCase)
            filePath = fullfile(testCase.ProjectRoot, 'data_for_examples', ...
                'xlsx_templates', 'transducer_sf.xlsx');

            actual = read_transducer_sf_from_xls(filePath);

            testCase.verifyEqual(actual.expSign, 1);
            testCase.verifyGreaterThan(actual.frequency, 0);
            testCase.verifySize(actual.xGrid, size(actual.complexVelocityAmplitude));
            testCase.verifyTrue(all(isfinite(actual.complexVelocityAmplitude(:))));
        end

        function testHologramSpreadsheetTemplateLoads(testCase)
            filePath = fullfile(testCase.ProjectRoot, 'data_for_examples', ...
                'xlsx_templates', 'hologram_sf.xlsx');

            [geometry, hologram, medium] = read_hologram_sf_from_xls(filePath);

            testCase.verifyGreaterThan(geometry.apertureMin, 0);
            testCase.verifySize(hologram.xGrid, size(hologram.complexPressureAmplitude));
            testCase.verifyGreaterThan(medium.soundSpeed, 0);
            testCase.verifyGreaterThan(medium.density, 0);
        end

        function testExplicitCpuArchitecture(testCase)
            serviceParameters = struct('cpuArchitecture', 'sse2');

            [actual, isExplicit] = get_xddx_cpu_architecture(serviceParameters);

            testCase.verifyEqual(actual, 'sse2');
            testCase.verifyTrue(isExplicit);
        end

        function testExplicitCudaVersion(testCase)
            testCase.assumeFalse(ismac, 'CUDA is intentionally unsupported on macOS.');
            serviceParameters = struct('cudaVersion', 11);

            [actual, isExplicit] = get_xddx_cuda_version(serviceParameters);

            testCase.verifyEqual(actual, 'cuda11');
            testCase.verifyTrue(isExplicit);
        end

        function testSingleFrequencyCpuDockerImage(testCase)
            serviceParameters = struct('cpuArchitecture', 'avx2');

            actual = get_xddx_docker_image('cpu', serviceParameters, false);

            testCase.verifyEqual(actual, 'xddx-sf-cpu-avx2');
        end

        function testTransientCudaDockerImage(testCase)
            testCase.assumeFalse(ismac, 'CUDA is intentionally unsupported on macOS.');
            serviceParameters = struct('cudaVersion', 'cuda12');

            actual = get_xddx_docker_image('cuda', serviceParameters, true);

            testCase.verifyEqual(actual, 'xddx-transient-cuda12');
        end

        function testExplicitDeviceResolutionDoesNotUseDetection(testCase)
            [device, selection] = resolve_auto_simulation_device('CPU');

            testCase.verifyEqual(device, 'cpu');
            testCase.verifyFalse(selection.usesAuto);
            testCase.verifyEqual(selection.selectedDevice, 'cpu');
        end

        function testDockerConfigurationIsComplete(testCase)
            actual = xddx_docker_config();

            testCase.verifyNotEmpty(actual.dockerUsername);
            testCase.verifyGreaterThan(actual.imageUpdatePeriodDays, 0);
        end

        function testCanonicalWindowsDockerVariableEnablesDocker(testCase)
            testCase.assumeTrue(ispc, 'The override only applies on Windows.');
            oldCanonicalValue = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
            oldLegacyValue = getenv('XDDX_USE_DOCKER');
            testCase.addTeardown(@() setenv( ...
                'XDDX_USE_DOCKER_ON_WINDOWS', oldCanonicalValue));
            testCase.addTeardown(@() setenv( ...
                'XDDX_USE_DOCKER', oldLegacyValue));
            setenv('XDDX_USE_DOCKER_ON_WINDOWS', 'YeS');
            setenv('XDDX_USE_DOCKER', '0');

            actual = use_xddx_docker_on_windows();

            testCase.verifyTrue(actual);
        end

        function testCanonicalWindowsDockerVariableTakesPrecedence(testCase)
            testCase.assumeTrue(ispc, 'The override only applies on Windows.');
            oldCanonicalValue = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
            oldLegacyValue = getenv('XDDX_USE_DOCKER');
            testCase.addTeardown(@() setenv( ...
                'XDDX_USE_DOCKER_ON_WINDOWS', oldCanonicalValue));
            testCase.addTeardown(@() setenv( ...
                'XDDX_USE_DOCKER', oldLegacyValue));
            setenv('XDDX_USE_DOCKER_ON_WINDOWS', '0');
            setenv('XDDX_USE_DOCKER', '1');

            actual = use_xddx_docker_on_windows();

            testCase.verifyFalse(actual);
        end

        function testLegacyWindowsDockerVariableRemainsSupported(testCase)
            testCase.assumeTrue(ispc, 'The override only applies on Windows.');
            oldCanonicalValue = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
            oldLegacyValue = getenv('XDDX_USE_DOCKER');
            testCase.addTeardown(@() setenv( ...
                'XDDX_USE_DOCKER_ON_WINDOWS', oldCanonicalValue));
            testCase.addTeardown(@() setenv( ...
                'XDDX_USE_DOCKER', oldLegacyValue));
            setenv('XDDX_USE_DOCKER_ON_WINDOWS', '');
            setenv('XDDX_USE_DOCKER', 'true');

            actual = use_xddx_docker_on_windows();

            testCase.verifyTrue(actual);
        end
        function testDockerDaemonFailureHasActionableHint(testCase)
            output = 'Cannot connect to the Docker daemon. Is the docker daemon running?';

            actual = get_docker_run_failure_hint(output);

            testCase.verifySubstring(actual, 'Start Docker');
        end

        function testDockerImageFailureHasPullHint(testCase)
            output = 'Unable to find image locally: no such image';

            actual = get_docker_run_failure_hint(output);

            testCase.verifySubstring(actual, 'docker pull');
        end
    end
end
