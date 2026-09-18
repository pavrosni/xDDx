classdef XddxNumericalUtilitiesTest < matlab.unittest.TestCase
    %XDDXNUMERICALUTILITIESTEST Fast tests for reusable numerical helpers.

    methods (TestClassSetup)
        function addLibraryPath(testCase)
            projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            libraryDir = fullfile(projectRoot, 'xDDx_lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(libraryDir));
        end
    end

    methods (Test)
        function testShiftPiWrapsAtExpectedBoundaries(testCase)
            inputAngles = [-4*pi, -3*pi, -pi, 0, pi, 3*pi, 4*pi];
            expected = [0, -pi, -pi, 0, pi, pi, 0];

            actual = shift_pi(inputAngles);

            testCase.verifyEqual(actual, expected, 'AbsTol', 10*eps(pi));
        end

        function testShiftTwoPiPreservesPositiveMultiples(testCase)
            inputAngles = [-4*pi, -2*pi, -pi, 0, pi, 2*pi, 4*pi];
            expected = [0, 0, pi, 0, pi, 2*pi, 2*pi];

            actual = shift_2pi(inputAngles);

            testCase.verifyEqual(actual, expected, 'AbsTol', 10*eps(pi));
        end

        function testCenteredGridSupportsOddDimensions(testCase)
            [xGrid, yGrid] = build_flat_grid_centered(5, 0.5, 3, 2);

            testCase.verifyEqual(xGrid(1, :), [-1, -0.5, 0, 0.5, 1], ...
                'AbsTol', eps);
            testCase.verifyEqual(yGrid(:, 1), [-2; 0; 2], 'AbsTol', eps);
            testCase.verifyTrue(is_meshgrid(xGrid, yGrid));
        end

        function testCenteredGridSupportsEvenDimensions(testCase)
            [xGrid, yGrid] = build_flat_grid_centered(4, 0.5, 2, 2);

            testCase.verifyEqual(xGrid(1, :), [-1, -0.5, 0, 0.5], ...
                'AbsTol', eps);
            testCase.verifyEqual(yGrid(:, 1), [-2; 0], 'AbsTol', eps);
        end

        function testMeshgridRejectsNonCartesianCoordinates(testCase)
            xGrid = [0, 1; 0.25, 1];
            yGrid = [0, 0; 1, 1];

            actual = is_meshgrid(xGrid, yGrid);

            testCase.verifyFalse(actual);
        end

        function testWaveNumberArrayEvenLength(testCase)
            [waveNumbers, step] = make_k_array(4, 0.5);

            testCase.verifyEqual(step, pi, 'AbsTol', eps(pi));
            testCase.verifyEqual(waveNumbers, [-2, -1, 0, 1]*pi, ...
                'AbsTol', 10*eps(pi));
        end

        function testWaveNumberArrayOddLength(testCase)
            [waveNumbers, step] = make_k_array(5, 0.4);

            testCase.verifyEqual(step, pi, 'AbsTol', eps(pi));
            testCase.verifyEqual(waveNumbers, [-2, -1, 0, 1, 2]*pi, ...
                'AbsTol', 10*eps(pi));
        end

        function testFlatPressureVelocityForPlaneWave(testCase)
            medium = struct('soundSpeed', 1500, 'density', 1000);
            pressure = 2.5*ones(8, 6);
            expected = pressure/(medium.density*medium.soundSpeed);

            actual = pressure_to_velocity_flat_surface(pressure, 1e-3, ...
                1e-3, 1e6, medium);

            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-18);
        end

        function testHologramPowerIsFiniteAndPositive(testCase)
            medium = struct('soundSpeed', 1500, 'density', 1000);
            pressure = ones(16, 12);

            actual = holo_pressure_to_power(pressure, 0.25e-3, ...
                0.25e-3, 1e6, medium);

            testCase.verifyGreaterThan(actual, 0);
            testCase.verifyTrue(isfinite(actual));
        end

        function testReshapeTransducerConvertsVectorFields(testCase)
            input = struct('xGrid', 1:3, 'yGrid', 4:6, 'zGrid', 7:9, ...
                'complexPressureAmplitude', [1, 2, 3]);

            actual = reshape_transducer_dim(input);

            testCase.verifySize(actual.xGrid, [1, 1, 3]);
            testCase.verifySize(actual.yGrid, [1, 1, 3]);
            testCase.verifySize(actual.zGrid, [1, 1, 3]);
            testCase.verifySize(actual.complexPressureAmplitude, [1, 1, 3]);
        end

        function testStructureConsistencyAllowsTransientFrequencyAxis(testCase)
            input = struct('xGrid', zeros(2, 3), 'yGrid', zeros(2, 3), ...
                'zGrid', zeros(2, 3), 'input', zeros(2, 3, 4), ...
                'dx', 1, 'dy', 1);

            actual = check_consistency_of_struct(input, true);

            testCase.verifyTrue(actual);
        end

        function testStructureConsistencyRejectsMismatchedGrids(testCase)
            input = struct('xGrid', zeros(2, 3), 'yGrid', zeros(3, 2), ...
                'zGrid', zeros(2, 3));

            actual = check_consistency_of_struct(input, false);

            testCase.verifyFalse(actual);
        end

        function testValidForwardPlaneInputIsAccepted(testCase)
            source = struct('xGrid', [0; 1e-3], 'yGrid', [0; 0], ...
                'zGrid', [0; 0], 'dx', 1e-3, 'dy', 1e-3, ...
                'input', [1; 2]);
            field = struct('xGrid', 0, 'yGrid', 0, 'zGrid', 1e-2);
            messages = load_error_messages();

            [input, elementArea, isSource, isField] = check_input_errors( ...
                'cpu', false, 4, 1:6, messages, source, field, ...
                eps('single'), 1);

            testCase.verifyEqual(input, [1; 2], 'AbsTol', eps);
            testCase.verifyEqual(elementArea, 1e-6, 'AbsTol', eps(1e-6));
            testCase.verifyTrue(isSource);
            testCase.verifyFalse(isField);
        end

        function testFlatTransducerDisablesAlignment(testCase)
            [doAlignment, warningMessage] = check_flat_transducer_issue(false, true);

            testCase.verifyFalse(doAlignment);
            testCase.verifyNotEmpty(warningMessage);
        end

        function testMatlabVideoOptionsRemainEnabled(testCase)
            [saveSpectrum, saveSignal] = check_octave_based_video_issues( ...
                false, true, true);

            testCase.verifyTrue(saveSpectrum);
            testCase.verifyTrue(saveSignal);
        end
    end
end
