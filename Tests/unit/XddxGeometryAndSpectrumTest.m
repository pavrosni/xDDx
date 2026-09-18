classdef XddxGeometryAndSpectrumTest < matlab.unittest.TestCase
    %XDDXGEOMETRYANDSPECTRUMTEST Geometry, FFT, and plotting contracts.

    methods (TestClassSetup)
        function addLibraryPath(testCase)
            projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            libraryDir = fullfile(projectRoot, 'xDDx_lib');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(libraryDir));
        end
    end

    methods (TestMethodSetup)
        function hideFigures(testCase)
            oldVisibility = get(groot, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot, 'DefaultFigureVisible', oldVisibility));
            testCase.addTeardown(@() close('all', 'force'));
            set(groot, 'DefaultFigureVisible', 'off');
        end
    end

    methods (Test)
        function testPlaneFitUsesTolerance(testCase)
            zGrid = [1, 1 + 1e-9];

            withinTolerance = check_plane_fit(zGrid, 1e-8);
            outsideTolerance = check_plane_fit(zGrid, 1e-10);

            testCase.verifyTrue(withinTolerance);
            testCase.verifyFalse(outsideTolerance);
        end

        function testSphereFitRecognizesSphericalCap(testCase)
            radius = 0.08;
            [xGrid, yGrid] = meshgrid([-0.01, 0, 0.01]);
            zGrid = radius - sqrt(radius^2 - xGrid.^2 - yGrid.^2);

            actual = check_sphere_fit(xGrid, yGrid, zGrid, 1e-12, radius);

            testCase.verifyTrue(actual);
        end

        function testTransducerZDefaultsToFlatSurface(testCase)
            input = struct('xGrid', [-1, 0, 1], 'yGrid', [0, 0, 0]);

            actual = check_transducer_z(input);

            testCase.verifyEqual(actual.zGrid, zeros(1, 3), 'AbsTol', eps);
        end

        function testTransducerZDefaultsToSphericalCap(testCase)
            radius = 2;
            input = struct('xGrid', [-1, 0, 1], 'yGrid', [0, 0, 0], ...
                'radiusOfCurvature', radius);
            expected = radius - sqrt(radius^2 - input.xGrid.^2);

            actual = check_transducer_z(input);

            testCase.verifyEqual(actual.zGrid, expected, 'AbsTol', 10*eps(radius));
        end

        function testRotatedCoordinatesForAxialDirection(testCase)
            xGrid = [-0.01, 0.01];
            yGrid = [0, 0];

            [xRotated, yRotated, zRotated, zPosition] = ...
                calculate_rotated_coordinates(0, 0, 0.1, [0; 0; 1], ...
                xGrid, yGrid, 0.1, 0.01);

            testCase.verifyEqual(xRotated, xGrid, 'AbsTol', eps);
            testCase.verifyEqual(yRotated, yGrid, 'AbsTol', eps);
            testCase.verifyEqual(zRotated, 0.02*ones(size(xGrid)), ...
                'AbsTol', 10*eps(0.02));
            testCase.verifyEqual(zPosition, 0.02, 'AbsTol', 10*eps(0.02));
        end

        function testSingleSidedSpectrumLocatesTone(testCase)
            sampleRate = 8e3;
            sampleCount = 64;
            toneFrequency = 1e3;
            time = (0:(sampleCount - 1))/sampleRate;
            waveform = reshape(cos(2*pi*toneFrequency*time + pi/5), 1, 1, []);

            [frequencies, spectrum] = get_single_sided_spectrum(time, waveform);
            [~, peakIndex] = max(abs(spectrum), [], 3);

            testCase.verifyEqual(frequencies(peakIndex), toneFrequency, 'AbsTol', eps(toneFrequency));
            testCase.verifyEqual(abs(spectrum(peakIndex)), 1, 'AbsTol', 1e-12);
        end

        function testExtractSingleFrequencyMatchesSpectrum(testCase)
            sampleRate = 8e3;
            sampleCount = 64;
            toneFrequency = 1e3;
            time = (0:(sampleCount - 1))/sampleRate;
            waveform = reshape(cos(2*pi*toneFrequency*time + pi/5), 1, 1, []);
            transient = struct('time', time, 'pressureWaveforms', waveform, ...
                'xGrid', 0, 'yGrid', 0, 'zPosition', 0.02, ...
                'dx', 1e-3, 'dy', 1e-3);
            [frequencies, spectrum] = get_single_sided_spectrum(time, waveform);
            [~, expectedIndex] = min(abs(frequencies - toneFrequency));

            [hologram, actualFrequency] = extract_sf(transient, 1, toneFrequency);

            testCase.verifyEqual(actualFrequency, toneFrequency, 'AbsTol', eps(toneFrequency));
            testCase.verifyEqual(hologram.complexPressureAmplitude, ...
                spectrum(:, :, expectedIndex), 'AbsTol', 1e-12);
        end

        function testSignificantFrequencyBounds(testCase)
            spectrum = reshape([0, 0.1, 1, 0.5, 0.01], 1, 1, []);

            lowest = find_lowest_significant_freq(spectrum, 0.04);
            highest = find_highest_significant_freq(100, spectrum, 0.04, 20, 1e6);
            estimated = estimate_highest_significant_freq(100, spectrum, 0.04, 20, 1e6);

            testCase.verifyEqual(lowest, 2);
            testCase.verifyEqual(highest, 4);
            testCase.verifyEqual(estimated, 4);
        end

        function testWaveformCenteringMovesPeakAndReturnsBounds(testCase)
            waveform = zeros(2, 2, 12);
            waveform(:, :, 3) = 1;

            [centered, firstIndex, lastIndex] = circshift_waveform( ...
                waveform, 0.1, false);
            peakByTime = squeeze(max(max(abs(centered), [], 1), [], 2));
            [~, peakIndex] = max(peakByTime);

            testCase.verifyEqual(peakIndex, 5);
            testCase.verifyGreaterThanOrEqual(firstIndex, 1);
            testCase.verifyLessThanOrEqual(lastIndex, size(waveform, 3));
        end

        function testHologramRotationVectorRecoversAxialLine(testCase)
            [xGrid, yGrid, zGrid] = meshgrid([-1, 0, 1], [-1, 0, 1], 1:4);
            pressure = zeros(size(xGrid));
            pressure(2, 2, :) = 1:4;

            [xMax, yMax, zMax, direction, line] = ...
                hologram_rotation_vector(xGrid, yGrid, zGrid, pressure, 1:4);

            testCase.verifyEqual([xMax, yMax, zMax], [0, 0, 4], 'AbsTol', eps);
            testCase.verifyEqual(direction, [0; 0; 1], 'AbsTol', 1e-12);
            testCase.verifyEqual(line.angleZ, 0, 'AbsTol', 1e-12);
            testCase.verifySize(line.xyzInputPoints, [3, 4]);
        end

        function testDisplayZeroDimensionalFieldReturnsInputs(testCase)
            outputText = evalc(['[xOut, yOut, zOut, pOut] = ' ...
                'disp_0d_field(1e-3, 2e-3, 3e-3, 4 + 3i);']);

            testCase.verifyEqual(xOut, 1e-3, 'AbsTol', eps);
            testCase.verifyEqual(yOut, 2e-3, 'AbsTol', eps);
            testCase.verifyEqual(zOut, 3e-3, 'AbsTol', eps);
            testCase.verifyEqual(pOut, 4 + 3i, 'AbsTol', eps);
            testCase.verifySubstring(outputText, '5 Pa');
        end

        function testPlotOneDimensionalFieldReturnsInputData(testCase)
            xGrid = reshape(0:3, 1, 1, []);
            yGrid = zeros(size(xGrid));
            zGrid = zeros(size(xGrid));
            pressure = reshape(1:4, 1, 1, []);

            [xOut, yOut, zOut, pOut] = plot_1d_field( ...
                xGrid, yGrid, zGrid, pressure);

            testCase.verifyEqual(xOut, squeeze(xGrid), 'AbsTol', eps);
            testCase.verifyEqual(yOut, squeeze(yGrid), 'AbsTol', eps);
            testCase.verifyEqual(zOut, squeeze(zGrid), 'AbsTol', eps);
            testCase.verifyEqual(pOut, squeeze(pressure), 'AbsTol', eps);
        end

        function testPlotTwoDimensionalFieldReturnsInputData(testCase)
            [xGrid, yGrid] = meshgrid(0:2, 0:1);
            zGrid = zeros(size(xGrid));
            pressure = xGrid + 2*yGrid;

            [xOut, yOut, zOut, pOut] = plot_2d_field( ...
                xGrid, yGrid, zGrid, pressure);

            testCase.verifyEqual(xOut, xGrid, 'AbsTol', eps);
            testCase.verifyEqual(yOut, yGrid, 'AbsTol', eps);
            testCase.verifyEqual(zOut, zGrid, 'AbsTol', eps);
            testCase.verifyEqual(pOut, pressure, 'AbsTol', eps);
        end
    end
end
