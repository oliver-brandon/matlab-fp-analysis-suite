choice = input(:,1);
reward = input(:,2);

for i = 1:size(choice,1)
    if choice(i) == 1
        cueCorrect(i,:) = sigs(i,:);
    elseif choice(i) == 0
        cueIncorrect(i,:) = sigs(i,:);
    end
end

meanCueCorrect = mean(cueCorrect,1);
meanCueIncorrect = mean(cueIncorrect,1);

numCueCorDiv = round(size(cueCorrect,1)/2);
numCueIncDiv = round(size(cueIncorrect,1)/2);

meanCueCorFirst = mean(cueCorrect(1:numCueCorDiv,:),1);
meanCueCorSecond = mean(cueCorrect(numCueCorDiv+1:end,:),1);

meanCueIncFirst = mean(cueIncorrect(1:numCueIncDiv,:),1);
meanCueIncSecond = mean(cueIncorrect(numCueIncDiv+1:end,:),1);

input = [];
sigs = [];
