clearvars -except input

path = ...
'/Users/brandon/personal-drive/optomouse-prime/opto-reversal/cue-stim/mats-cuestim-VTA/Ipsilateral/';
file = '1247M_Rev3CueStim_Chrim.mat';
filename = fullfile(path, file);
S = load(filename);
data = S.data;

leverOnset  = data.epocs.levers.onset;
rewardOnset = input;   % expected size: [numel(leverOnset) x 2]
cRew = [];
cNoRew = [];
iRew = [];
iNoRew = [];

n = numel(leverOnset);
if size(rewardOnset,1) ~= n || size(rewardOnset,2) < 2
    error('rewardOnset must be an N x 2 matrix where N = numel(leverOnset).');
end

for i = 1:n
    t = leverOnset(i);   % time of the lever press for this trial

    isCorrect = rewardOnset(i,1) == 1;
    isReward  = rewardOnset(i,2) == 1;

    if isCorrect && isReward
        cRew = [cRew; t];
    elseif isCorrect && ~isReward
        cNoRew = [cNoRew; t];
    elseif ~isCorrect && isReward
        iRew = [iRew; t];
    else
        iNoRew = [iNoRew; t];
    end
end

if ~isempty(cRew)
    data.epocs.cRewA.onset  = cRew;
    data.epocs.cRewA.offset = cRew + 1;
    data.epocs.cRewA.name   = 'cRewA';
    data.epocs.cRewA.data   = ones(numel(cRew), 1);
end
if ~isempty(cNoRew)
    data.epocs.cNoRewA.onset  = cNoRew;
    data.epocs.cNoRewA.offset = cNoRew + 1;
    data.epocs.cNoRewA.name   = 'cNoRewA';
    data.epocs.cNoRewA.data   = ones(numel(cNoRew), 1);
end
if ~isempty(iRew)
    data.epocs.iRewA.onset  = iRew;
    data.epocs.iRewA.offset = iRew + 1;
    data.epocs.iRewA.name   = 'iRewA';
    data.epocs.iRewA.data   = ones(numel(iRew), 1);
end
if ~isempty(iNoRew)
    data.epocs.iNoRewA.onset  = iNoRew;
    data.epocs.iNoRewA.offset = iNoRew + 1;
    data.epocs.iNoRewA.name   = 'iNoRewA';
    data.epocs.iNoRewA.data   = ones(numel(iNoRew), 1);
end

data.epocs.events.data = rewardOnset;

save(filename,'data');