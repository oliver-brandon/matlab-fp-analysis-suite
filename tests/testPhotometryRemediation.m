function tests = testPhotometryRemediation
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repoRoot));
testCase.TestData.repoRoot = repoRoot;
end

function testCorrectSignalDirection(testCase)
t = (1:1000)';
F405 = 1 + 0.001 * t + 0.02 * sin(t / 25);
neural = 0.05 * exp(-((t - 500) / 60).^2);
F465 = 2 * F405 + 0.25 + neural;

corrected = correctSignal(F405, F465);

verifyLessThan(testCase, abs(mean(corrected(1:100))), 1e-10);
verifyGreaterThan(testCase, max(corrected), 0.04);
end

function testRewardWindowClassification(testCase)
lever = [10; 20; 30];
pellet = [10.25; 30.30];

[isRewarded, rewarded, unrewarded] = classifyRewardedLevers(lever, pellet, 2);

verifyEqual(testCase, isRewarded, [true; false; true]);
verifyEqual(testCase, rewarded, [10; 30]);
verifyEqual(testCase, unrewarded, 20);
end

function testPhotometryChannelConfig(testCase)
cfgA = photometryChannelConfig('A');
cfgC = photometryChannelConfig(2);

verifyEqual(testCase, cfgA.signal, 'x465A');
verifyEqual(testCase, cfgC.correct, 'CL2_');
end

function testValidatePhotometryData(testCase)
data.streams.x405A.data = (1:100)';
data.streams.x405A.fs = 10;
data.epocs.St1_.onset = [1; 2];

[ok, problems] = validatePhotometryData(data, {'x405A'}, {'St1_'}, 5);
verifyTrue(testCase, ok);
verifyEmpty(testCase, problems);

[ok, problems] = validatePhotometryData(data, {'x465A'}, {'St1_'}, 5);
verifyFalse(testCase, ok);
verifyNotEmpty(testCase, problems);
end

function testPrlEpocsRewardWindow(testCase)
data = makeSyntheticPrlData();

data = prl_epocs(data);

verifyEqual(testCase, data.epocs.cRewA.onset, 11);
verifyEqual(testCase, data.epocs.cNoRewA.onset, 21);
verifyEqual(testCase, data.epocs.iRewA.onset, 31);
verifyEqual(testCase, data.epocs.iNoRewA.onset, 41);
verifyEqual(testCase, data.analysis.pipelineVersion, '2026-05-corrected');
end

function testPrlDfEpocsUsesChannelTwoTtls(testCase)
data = makeSyntheticPrlData();

data = prl_df_epocs(data, 2);

verifyEqual(testCase, data.epocs.cRewC.onset, 12);
verifyEqual(testCase, data.epocs.cNoRewC.onset, 22);
verifyEqual(testCase, data.epocs.iRewC.onset, 32);
verifyEqual(testCase, data.epocs.iNoRewC.onset, 42);
verifyEqual(testCase, data.epocs.levers.onset, [12; 22; 32; 42]);
end

function testErrorProbLoseShiftOffset(testCase)
data.epocs = struct();
sessionIdentifiers = [1 0; 2 2; 3 0; 4 3];

[data, errorProbLeverTS] = errorProbExtract(data, sessionIdentifiers, 4, 1);

verifyEqual(testCase, errorProbLeverTS, [2 4]);
verifyEqual(testCase, data.epocs.los_SHIFT_LEV.onset, 4);
verifyEqual(testCase, data.epocs.los_SHIFT_LEV.offset, 5);
end

function testProcessChunksIteratesRows(testCase)
arrValues = [
    0 0 5 0 0;
    0 0 0 6 0
];
arrIndexes = [
    1 2 3 4 5;
    6 7 8 9 10
];

[~, allIndices] = processChunks(arrValues, arrIndexes, 100, 1, 1, false);

verifyEqual(testCase, numel(allIndices), 2);
verifyEqual(testCase, allIndices{1}, 3);
verifyEqual(testCase, allIndices{2}, 4);
end

function testEpocExtractDoesNotShiftBaseline(testCase)
fs = 10;
sessionTime = 0:1/fs:10;
sessionSignal = 1 + 0.01 * sessionTime;
sessionSignal(sessionTime >= 5 & sessionTime <= 6) = 2;
TTLarray = 5;
ts1 = -2 + (1:70) / fs;

epocStream = epocExtract(sessionSignal, sessionTime, TTLarray, 2, 5, [-2 -1], [0 1], 70, ts1);
baselineIdx = ts1 >= -2 & ts1 <= -1;

verifyLessThan(testCase, abs(mean(epocStream(1,baselineIdx), 'omitnan')), 1e-10);
end

function data = makeSyntheticPrlData()
fs = 10;
sig = (1:1000)' / 1000;
data.streams.x405A.data = sig;
data.streams.x405A.fs = fs;
data.streams.x465A.data = sig;
data.streams.x465A.fs = fs;
data.streams.x405C.data = sig;
data.streams.x405C.fs = fs;
data.streams.x465C.data = sig;
data.streams.x465C.fs = fs;

data.epocs.St1_.onset = [10; 20; 30; 40];
data.epocs.CL1_.onset = [11; 21];
data.epocs.IL1_.onset = [31; 41];
data.epocs.Pe1_.onset = [11.25; 31.5];

data.epocs.St2_.onset = [10; 20; 30; 40];
data.epocs.CL2_.onset = [12; 22];
data.epocs.IL2_.onset = [32; 42];
data.epocs.Pe2_.onset = [12.25; 32.5];
end
