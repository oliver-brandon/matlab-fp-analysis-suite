function [data] = prl_epocs(data)

defaults = analysisDefaults();

if isfield(data.streams, 'x405A')
    cfg = photometryChannelConfig('A');
    emptyValue = 0;
elseif isfield(data.streams, 'x405C')
    cfg = photometryChannelConfig('C');
    emptyValue = 0;
else
    error('prl_epocs:MissingStreams', 'Cannot find x405A or x405C streams.');
end

data = addPrlEpocsForChannel(data, cfg, defaults.rewardWindowSec, emptyValue);
data.analysis.pipelineVersion = defaults.pipelineVersion;

end

function data = addPrlEpocsForChannel(data, cfg, rewardWindowSec, emptyValue)
requiredEpocs = {cfg.cue, cfg.pellet};
[ok, problems] = validatePhotometryData(data, {}, requiredEpocs, 0);
if ~ok
    error('prl_epocs:MissingEpocs', strjoin(problems, '; '));
end

cueTs = data.epocs.(cfg.cue).onset(:);
correctTs = getEpocOnsets(data, cfg.correct);
incorrectTs = getEpocOnsets(data, cfg.incorrect);
pelletTs = data.epocs.(cfg.pellet).onset(:);

[~, correctRewarded, correctNoReward] = classifyRewardedLevers(correctTs, pelletTs, rewardWindowSec);
[~, incorrectRewarded, incorrectNoReward] = classifyRewardedLevers(incorrectTs, pelletTs, rewardWindowSec);

correctRewardedOut = withEmptyValue(correctRewarded, emptyValue);
correctNoRewardOut = withEmptyValue(correctNoReward, emptyValue);
incorrectRewardedOut = withEmptyValue(incorrectRewarded, emptyValue);
incorrectNoRewardOut = withEmptyValue(incorrectNoReward, emptyValue);

data.epocs.(cfg.cRew) = makeEpoc(cfg.cRew, correctRewardedOut, 1);
data.epocs.(cfg.cNoRew) = makeEpoc(cfg.cNoRew, correctNoRewardOut, 2);
data.epocs.(cfg.iRew) = makeEpoc(cfg.iRew, incorrectRewardedOut, 3);
data.epocs.(cfg.iNoRew) = makeEpoc(cfg.iNoRew, incorrectNoRewardOut, 4);

correctLevers = sort([correctRewarded; correctNoReward], 'ascend');
incorrectLevers = sort([incorrectRewarded; incorrectNoReward], 'ascend');
cueCorrect = precedingCues(cueTs, correctLevers);
cueIncorrect = precedingCues(cueTs, incorrectLevers);

data.epocs.(cfg.cueCor) = makeEpoc(cfg.cueCor, cueCorrect, 5);
data.epocs.(cfg.cueInc) = makeEpoc(cfg.cueInc, cueIncorrect, 6);
end

function onsets = getEpocOnsets(data, epocName)
if isfield(data.epocs, epocName)
    onsets = data.epocs.(epocName).onset(:);
else
    onsets = [];
end
end

function values = withEmptyValue(values, emptyValue)
if isempty(values) && ~isempty(emptyValue)
    values = emptyValue;
end
values = values(:);
end

function epoc = makeEpoc(name, onsets, dataValue)
onsets = onsets(:);
epoc.name = name;
epoc.onset = onsets;
epoc.offset = onsets + 1;
epoc.data = ones(numel(onsets), 1) * dataValue;
end

function cues = precedingCues(cueTs, leverTs)
cues = [];
for i = 1:numel(leverTs)
    idx = find(cueTs < leverTs(i), 1, 'last');
    if ~isempty(idx)
        cues(end+1,1) = cueTs(idx); %#ok<AGROW>
    end
end
end
