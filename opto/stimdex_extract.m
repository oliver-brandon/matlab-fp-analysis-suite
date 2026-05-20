clearvars -except data
% if ~isfield(data.epocs, 'CL1_')
%     leverTS = data.epocs.IL1_.onset;
% elseif ~isfield(data.epocs, 'IL1_')
%     leverTS = data.epocs.CL1_.onset;
% elseif isfield(data.epocs, 'CL1_') && isfield(data.epocs, 'IL1_')
%     correct = data.epocs.CL1_.onset;
%     incorrect = data.epocs.IL1_.onset;
%     leverTS = sort([correct;incorrect]);
% end

leverTS = data.epocs.levers.onset;
stims = data.epocs.P1__.onset;
if stims > 10
    stims(11:end) = [];
end

for i = 1:height(stims)
    stimdex(i,1) = find(leverTS > stims(i,1),1);
end
