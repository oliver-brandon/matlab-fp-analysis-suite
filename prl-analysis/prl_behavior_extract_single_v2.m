clearvars -except choice outcome input;
% Create new variable called input first
% Right-click the workspace and click new variable OR
% in the command window, paste the following: input = [] and hit enter
choice = input(:,1);
outcome = input(:,2);
trials = length(choice);
per = 0;
reg = 0;
wS = 0;
lS = 0;
wS_inc = 0;
lS_inc = 0;
win_tot = 0;
lose_tot = 0;
win_tot_incorrect = 0;
lose_tot_incorrect = 0;
percent_wS = 0;
percent_lS = 0;

% Additional strategies (requested)
wShift = 0;
lStay = 0;
wShift_inc = 0;
lStay_inc = 0;
win_tot_next = 0;
lose_tot_next = 0;
win_tot_next_incorrect = 0;
lose_tot_next_incorrect = 0;
percent_wShift = 0;
percent_lStay = 0;
percent_wShift_inc = 0;
percent_lStay_inc = 0;

% Percent correct
percent_cor = sum(choice)/size(choice,1);
% perseveration
for i = 1:trials-2
    if choice(i) == 0
        per = per+1;
    end
    %if there are three 1s in a row (i+2) continue
    if choice(i) == 1 && choice(i+1) == 1 && choice(i+2) == 1
        tillthree = i+2;
        break
    else
        tillthree = 33;
        continue
    end
    if i == trials-2
        if choice(i) == 0
            per = per+1;
        end
        if choice(i+1) == 0
            per = per+1;
        end
        if choice(i+2) == 0
            per = per+1;
        end
    end
    
end
% per = per / sum(choice(1:end) == 0);
%regressive
%count the number of 0s after three 1s in a row. If there are no instances of three 1s in a row, set reg = 0


for i = 1:trials-2
    if choice(i) == 1 && choice(i+1) == 1 && choice(i+2) == 1
        reg = sum(choice(i+3:end) == 0);
        break
    end
    
        
end
% reg = reg / sum(choice(1:end) == 0);
%win-stay
for i = 1:trials-1
    if choice(i) == 1 && outcome(i) == 1
        win_tot = win_tot+1;
        if choice(i+1) == 1
            wS = wS+1;
        end
    end
    if i == trials-1
        if choice(i) == 1 && outcome(i) == 1
            win_tot = win_tot+1;
        end
    end
end
percent_wS = (wS/win_tot);
if isnan(percent_wS)
    percent_wS = 0;
else
    disp('')
end

%lose-shift
for i = 1:trials-1
    if choice(i) == 1 && outcome(i) == 0
        lose_tot = lose_tot+1;
        if choice(i+1) == 0
            lS = lS+1;
        end
    end
    if i == trials-1
        if choice(i) == 1 && outcome(i) == 0
            lose_tot = lose_tot+1;
        end
    end
end
percent_lS = (lS/lose_tot);
if isnan(percent_lS)
    percent_lS = 0;
else
    disp('')
end

% win-stay (incorrect lever)
for i = 1:trials-1
    if choice(i) == 0 && outcome(i) == 1
        win_tot_incorrect = win_tot_incorrect+1;
        if choice(i+1) == 0
            wS_inc = wS_inc+1;
        end
    end
    if i == trials-1
        if choice(i) == 0 && outcome(i) == 1
            win_tot_incorrect = win_tot_incorrect+1;
        end
    end
end
percent_wS_inc = (wS_inc/win_tot_incorrect);
if isnan(percent_wS_inc)
    percent_wS_inc = 0;
else
    disp('')
end
%lose-shift (incorrect lever)
for i = 1:trials-1
    if choice(i) == 0 && outcome(i) == 0
        lose_tot_incorrect = lose_tot_incorrect+1;
        if choice(i+1) == 1
            lS_inc = lS_inc+1;
        end
    end
    if i == trials-1
        if choice(i) == 0 && outcome(i) == 0
            lose_tot_incorrect = lose_tot_incorrect+1;
        end
    end
end
percent_lS_inc = (lS_inc/lose_tot_incorrect);
if isnan(percent_lS_inc)
    percent_lS_inc = 0;
else
    disp('')
end

% win-shift (correct lever): rewarded on correct lever, then switch to incorrect lever
for i = 1:trials-1
    if choice(i) == 1 && outcome(i) == 1
        win_tot_next = win_tot_next+1;
        if choice(i+1) == 0
            wShift = wShift+1;
        end
    end
end
percent_wShift = (wShift/win_tot_next);
if isnan(percent_wShift)
    percent_wShift = 0;
else
    disp('')
end

% lose-stay (correct lever): not rewarded on correct lever, then stay on correct lever
for i = 1:trials-1
    if choice(i) == 1 && outcome(i) == 0
        lose_tot_next = lose_tot_next+1;
        if choice(i+1) == 1
            lStay = lStay+1;
        end
    end
end
percent_lStay = (lStay/lose_tot_next);
if isnan(percent_lStay)
    percent_lStay = 0;
else
    disp('')
end

% win-shift (incorrect lever): rewarded on incorrect lever, then switch to correct lever
for i = 1:trials-1
    if choice(i) == 0 && outcome(i) == 1
        win_tot_next_incorrect = win_tot_next_incorrect+1;
        if choice(i+1) == 1
            wShift_inc = wShift_inc+1;
        end
    end
end
percent_wShift_inc = (wShift_inc/win_tot_next_incorrect);
if isnan(percent_wShift_inc)
    percent_wShift_inc = 0;
else
    disp('')
end

% lose-stay (incorrect lever): not rewarded on incorrect lever, then stay on incorrect lever
for i = 1:trials-1
    if choice(i) == 0 && outcome(i) == 0
        lose_tot_next_incorrect = lose_tot_next_incorrect+1;
        if choice(i+1) == 0
            lStay_inc = lStay_inc+1;
        end
    end
end
percent_lStay_inc = (lStay_inc/lose_tot_next_incorrect);
if isnan(percent_lStay_inc)
    percent_lStay_inc = 0;
else
    disp('')
end

fprintf('percent_cor = %d\nper = %d\nreg = %d\npercent_wS = %d\npercent_lS = %d\ntot trials = %d\ntill three = %d\npercentWS-inc = %d\npercentLS-inc = %d\n',...
    percent_cor,per,reg,percent_wS,percent_lS,trials,tillthree,percent_wS_inc,percent_lS_inc)
behaviorOutput = [percent_cor,per,reg,percent_wS,percent_lS,percent_wShift,percent_lStay,percent_wS_inc,percent_lS_inc,percent_wShift_inc,percent_lStay_inc];
outputTable = array2table(behaviorOutput,'VariableNames',{'percent_cor','per','reg','wS','lS','wShift','lStay','wS_inc','lS_inc','wShift_inc','lStay_inc'});
disp(head(outputTable))


%% --- Trials to 3 consecutive correct choices (choice==1) ---
% In Brandon's dataset, choice: 1 = correct, 0 = incorrect
trialsTo3Correct = NaN;
for i = 1:(numel(choice)-2)
    if choice(i)==1 && choice(i+1)==1 && choice(i+2)==1
        trialsTo3Correct = i+2;
        break;
    end
end

% if exist('behaviorOutput','var')
%     behaviorOutput.trialsTo3Correct = trialsTo3Correct;
% end
% 
% if exist('outputTable','var')
%     outputTable.trialsTo3Correct = trialsTo3Correct;
% end

fprintf('Trials to 3 consecutive correct: %g\n', trialsTo3Correct);

