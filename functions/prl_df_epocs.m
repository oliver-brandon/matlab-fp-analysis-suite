function [data] = prl_df_epocs(data, TTLs)

defaults = analysisDefaults();

if TTLs == 1
    cfg = photometryChannelConfig('A');
elseif TTLs == 2
    cfg = photometryChannelConfig('C');
else
    error('prl_df_epocs:InvalidTTLs', 'TTLs must be 1 for A or 2 for C.');
end

data = addPrlEpocsForChannel(data, cfg, defaults.rewardWindowSec);
data.analysis.pipelineVersion = defaults.pipelineVersion;

end

function data = addPrlEpocsForChannel(data, cfg, rewardWindowSec)
requiredEpocs = {cfg.cue, cfg.pellet};
[ok, problems] = validatePhotometryData(data, {}, requiredEpocs, 0);
if ~ok
    error('prl_df_epocs:MissingEpocs', strjoin(problems, '; '));
end

cueTs = data.epocs.(cfg.cue).onset(:);
correctTs = getEpocOnsets(data, cfg.correct);
incorrectTs = getEpocOnsets(data, cfg.incorrect);
pelletTs = data.epocs.(cfg.pellet).onset(:);

[~, correctRewarded, correctNoReward] = classifyRewardedLevers(correctTs, pelletTs, rewardWindowSec);
[~, incorrectRewarded, incorrectNoReward] = classifyRewardedLevers(incorrectTs, pelletTs, rewardWindowSec);

data.epocs.(cfg.cRew) = makeEpoc(cfg.cRew, correctRewarded, 1);
data.epocs.(cfg.cNoRew) = makeEpoc(cfg.cNoRew, correctNoReward, 2);
data.epocs.(cfg.iRew) = makeEpoc(cfg.iRew, incorrectRewarded, 3);
data.epocs.(cfg.iNoRew) = makeEpoc(cfg.iNoRew, incorrectNoReward, 4);

correctLevers = sort([correctRewarded; correctNoReward], 'ascend');
incorrectLevers = sort([incorrectRewarded; incorrectNoReward], 'ascend');
cueCorrect = precedingCues(cueTs, correctLevers);
cueIncorrect = precedingCues(cueTs, incorrectLevers);

data.epocs.(cfg.cueCor) = makeEpoc(cfg.cueCor, cueCorrect, 5);
data.epocs.(cfg.cueInc) = makeEpoc(cfg.cueInc, cueIncorrect, 6);

levers = sort([correctTs; incorrectTs]);
data.epocs.levers = makeEpoc('levers', levers, 1);
end

function onsets = getEpocOnsets(data, epocName)
if isfield(data.epocs, epocName)
    onsets = data.epocs.(epocName).onset(:);
else
    onsets = [];
end
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
