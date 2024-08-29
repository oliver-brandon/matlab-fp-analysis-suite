clearvars -except choice outcome;

% choice = [0;0;0;1;1;1;1;1;1;0;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
% outcome = [0;0;0;1;0;1;1;1;0;1;1;1;1;1;0;0;1;1;1;1;1;0;1;1;1;0;1;1;1;1];

trials = length(choice);
per = 0;
reg = 0;
wS = 0;
lS = 0;
win_tot = 0;
lose_tot = 0;
percent_wS = 0;
percent_lS = 0;
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

%regressive
%count the number of 0s after three 1s in a row. If there are no instances of three 1s in a row, set reg = 0


for i = 1:trials-2
    if choice(i) == 1 && choice(i+1) == 1 && choice(i+2) == 1
        reg = sum(choice(i+3:end) == 0);
    end
        
end

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

fprintf('per = %d\nreg = %d\npercent_wS = %d\npercent_lS = %d\ntot trials = %d\ntill three = %d\n',...
    per,reg,percent_wS,percent_lS,trials,tillthree)