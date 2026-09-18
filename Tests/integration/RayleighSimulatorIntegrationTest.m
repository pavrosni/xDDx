classdef RayleighSimulatorIntegrationTest < matlab.unittest.TestCase
    %RAYLEIGHSIMULATORINTEGRATIONTEST Tiny native CPU end-to-end cases.
    %
    % These tests exercise MATLAB-to-binary serialization, native executable
    % selection, output loading, and cleanup without using a large example.

    methods (TestClassSetup)
        function addLibraryPath(testCase)
            projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            libraryDir = fullfile(projectRoot, 'xDDx_lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(libraryDir));
        end
    end

    methods (TestMethodSetup)
        function forceNativeExecution(testCase)
            oldValue = getenv('XDDX_USE_DOCKER_ON_WINDOWS');
            testCase.addTeardown(@() setenv('XDDX_USE_DOCKER_ON_WINDOWS', oldValue));
            testCase.addTeardown(@() cd(testCase.ProjectRoot));
            if ~ispc
                dockerCheck = check_docker_ready();
                testCase.assumeTrue(dockerCheck.ok, ...
                    'Native executables are Windows-only; Docker is required on Unix/macOS.');
            end
            setenv('XDDX_USE_DOCKER_ON_WINDOWS', '0');
        end
    end

    methods (Test)
        function testSingleFrequencyCpuProjection(testCase)
            source = RayleighSimulatorIntegrationTest.singleFrequencySource();
            field = RayleighSimulatorIntegrationTest.fieldPoints();
            medium = RayleighSimulatorIntegrationTest.medium();
            service = struct('cpuArchitecture', get_xddx_cpu_architecture());

            actual = rayleigh_simulator(1, 1e6, 4, 'cpu', false, ...
                source, field, medium, service);

            testCase.verifySize(actual, size(field.xGrid));
            testCase.verifyTrue(all(isfinite(actual(:))));
            testCase.verifyGreaterThan(max(abs(actual(:))), 0);
        end

        function testTransientCpuProjection(testCase)
            source = RayleighSimulatorIntegrationTest.transientSource();
            field = RayleighSimulatorIntegrationTest.fieldPoints();
            medium = RayleighSimulatorIntegrationTest.medium();
            service = struct('cpuArchitecture', get_xddx_cpu_architecture());

            actual = rayleigh_simulator(1, 0.5e6, 4, 'cpu', true, ...
                source, field, medium, service);

            testCase.verifySize(actual, [size(field.xGrid), 2]);
            testCase.verifyTrue(all(isfinite(actual(:))));
            testCase.verifyGreaterThan(max(abs(actual(:))), 0);
        end
    end

    methods (Access=private)
        function root = ProjectRoot(~)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        end
    end

    methods (Static, Access=private)
        function source = singleFrequencySource()
            [xGrid, yGrid] = build_flat_grid_centered(3, 1e-3, 3, 1e-3);
            source = struct('xGrid', xGrid, 'yGrid', yGrid, ...
                'zGrid', zeros(size(xGrid)), 'dx', 1e-3, 'dy', 1e-3, ...
                'input', ones(size(xGrid))/(1500*1000));
        end

        function source = transientSource()
            source = RayleighSimulatorIntegrationTest.singleFrequencySource();
            source.input = cat(3, source.input, 0.5*source.input);
        end

        function field = fieldPoints()
            field = struct('xGrid', reshape([0, 1e-3], 1, 1, []), ...
                'yGrid', zeros(1, 1, 2), ...
                'zGrid', reshape([0.02, 0.03], 1, 1, []));
        end

        function value = medium()
            value = struct('soundSpeed', 1500, 'density', 1000);
        end
    end
end
