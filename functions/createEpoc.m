function [data] = createEpoc(data, epocOnset, epocName)

data.epocs.(epocName).onset = epocOnset;
data.epocs.(epocName).offset = epocOnset + 1;
data.epocs.(epocName).name = epocName;
data.epocs.(epocName).data = ones(height(epocOnset),1);

end