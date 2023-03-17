function [LOAD_BAR] = basicwaitbar(numFiles,i)

while i/numFiles < 0.25
    LOAD_BAR = waitbar(0,'Please Wait...');
end
while i/numFiles >= 0.25
    waitbar(0.25,LOAD_BAR,'25% Complete...');
end
while i/numFiles >= 0.50
    waitbar(0.5,LOAD_BAR,'50% Complete...');
end
while i/numFiles >= 0.75
    waitbar(0.75,LOAD_BAR,'75% Complete...');
end
while i/numFiles == 1
    waitbar(1,LOAD_BAR,'100% Complete!');
    close(LOAD_BAR)
end

