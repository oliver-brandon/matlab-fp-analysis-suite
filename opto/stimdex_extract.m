clearvars -except data
if ~isfield(data.epocs, 'CL1_')
    leverTS = data.epocs.IL1_.onset;
elseif ~isfield(data.epocs, 'IL1_')
    leverTS = data.epocs.CL1_.onset;
elseif isfield(data.epocs, 'CL1_') && isfield(data.epocs, 'IL1_')
    correct = data.epocs.CL1_.onset;
    incorrect = data.epocs.IL1_.onset;
    leverTS = sort([correct;incorrect]);
end

stims = data.epocs.P1__.onset;

for i = 1:height(stims)
    stimdex(i,1) = find(leverTS > stims(i,1),1);
end
