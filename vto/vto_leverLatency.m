clear;
blockpath = '/Users/brandon/personal-drive/vto/tanks/1996_VTOLC';

data = TDTbin2mat(blockpath);

cueOFF = data.epocs.St1_.offset;
leverON = data.epocs.RL1_.onset;

if size(cueOFF,1) > 30
    cueOFF = cueOFF(1:30,:);
end

if size(leverON,1) > 30
    leverON = leverON(1:30,:);
end

% Calculate the duration between cueOFF and leverON events
leverLatency = leverON - cueOFF;
meanLeverLatency = mean(leverLatency);