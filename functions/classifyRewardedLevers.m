function [isRewarded, rewardedLevers, unrewardedLevers] = classifyRewardedLevers(leverTs, pelletTs, rewardWindowSec)
% classifyRewardedLevers Classify lever timestamps by following pellet TTLs.

if nargin < 3 || isempty(rewardWindowSec)
    defaults = analysisDefaults();
    rewardWindowSec = defaults.rewardWindowSec;
end

leverTs = normalizeTs(leverTs);
pelletTs = normalizeTs(pelletTs);
isRewarded = false(size(leverTs));

if isempty(leverTs)
    rewardedLevers = leverTs;
    unrewardedLevers = leverTs;
    return
end

for i = 1:numel(leverTs)
    isRewarded(i) = any(pelletTs >= leverTs(i) & pelletTs <= leverTs(i) + rewardWindowSec);
end

rewardedLevers = leverTs(isRewarded);
unrewardedLevers = leverTs(~isRewarded);

end

function ts = normalizeTs(ts)
if isempty(ts) || isequal(ts, 0)
    ts = [];
else
    ts = ts(:);
    ts = ts(~isnan(ts) & ts ~= 0);
end
end
