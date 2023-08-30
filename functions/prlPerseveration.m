function [totPersev] = prlPerseveration(session_identifiers)

leverArray = session_identifiers(2:2:end,2);
patterns = [1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2; 2 2 1; 2 2 2];
idx = [];
for i = 1:height(leverArray) - 2
    for j = 1:size(patterns,1)
        if isequal(leverArray(i:i+2), patterns(j,:)')
            idx = [idx; i+2];
            break
        end
    end
end

if isempty(idx)
    totPersev = 33;
else
    totPersev = idx;
end

