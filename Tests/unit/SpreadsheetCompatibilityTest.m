classdef SpreadsheetCompatibilityTest < matlab.unittest.TestCase
    %SPREADSHEETCOMPATIBILITYTEST Legacy spreadsheet import edge cases.

    properties (SetAccess=private)
        TemporaryFolder
    end

    methods (TestClassSetup)
        function addLibraryPath(testCase)
            projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'xDDx_lib')));
        end
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.TemporaryFolder = fixture.Folder;
        end
    end

    methods (Test)
        function testBlankRadiusAndZProduceFlatTransducer(testCase)
            filePath = testCase.createWorkbook(NaN, {0, 0.001; 0, 0.001}, {NaN});

            actual = read_transducer_sf_from_xls(filePath);

            testCase.verifyEmpty(actual.radiusOfCurvature);
            testCase.verifyEqual(actual.zGrid, zeros(2), 'AbsTol', eps);
            testCase.verifyEqual(actual.complexVelocityAmplitude, ...
                [1, 2; 3, 4].*exp(1i*[0, 0.1; 0.2, 0.3]), 'AbsTol', 1e-12);
        end

        function testBlankZProducesSphericalTransducer(testCase)
            filePath = testCase.createWorkbook(0.08, {0, 0.001; 0, 0.001}, {NaN});

            actual = read_transducer_sf_from_xls(filePath);

            expectedZ = 0.08 - sqrt(0.08^2 - actual.xGrid.^2 - actual.yGrid.^2);
            testCase.verifyEqual(actual.zGrid, expectedZ, 'AbsTol', 1e-12);
        end

        function testNumericTextCellsRemainSupported(testCase)
            filePath = testCase.createWorkbook(NaN, {'0', '0.001'; '0', '0.001'}, {NaN});

            actual = read_transducer_sf_from_xls(filePath);

            testCase.verifyEqual(actual.xGrid, [0, 0.001; 0, 0.001], 'AbsTol', eps);
        end

        function testInvalidTextBorderIsRejected(testCase)
            filePath = testCase.createWorkbook(NaN, ...
                {0, 0.001, 'invalid'; 0, 0.001, 'invalid'}, {NaN});

            testCase.verifyError(@() read_transducer_sf_from_xls(filePath), '');
        end

        function testMissingInteriorCellIsRejected(testCase)
            filePath = testCase.createWorkbook(NaN, ...
                {0, NaN, 0.002; 0, 0.001, 0.002}, {NaN});

            testCase.verifyError(@() read_transducer_sf_from_xls(filePath), '');
        end
    end

    methods (Access=private)
        function filePath = createWorkbook(testCase, radius, xCells, zCells)
            filePath = fullfile(testCase.TemporaryFolder, 'transducer.xlsx');
            sheetNames = {'Info', 'Scalar Parameters', 'xGrid (m)', ...
                'yGrid (m)', 'zGrid (m)', 'velocityAmplitude (m_s)', 'velocityPhase (rad)'};
            sheets = {{'Compatibility fixture'}, ...
                {'expSign', 1; 'frequency', 1e6; 'radiusOfCurvature', radius; ...
                 'dx', 0.001; 'dy', 0.001}, ...
                xCells, {0, 0; 0.001, 0.001}, zCells, ...
                {1, 2; 3, 4}, {0, 0.1; 0.2, 0.3}};
            for sheetIndex = 1:numel(sheets)
                writetable(cell2table(sheets{sheetIndex}), filePath, ...
                    'Sheet', sheetNames{sheetIndex}, 'WriteVariableNames', false);
            end
        end
    end
end
