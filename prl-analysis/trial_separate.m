clearvars -except input;
close all;

correct_firsthalf = [];
correct_secondhalf = [];
incorrect_firsthalf = [];
incorrect_secondhalf = [];

for i = 1:height(input(:,1))
    if input(i,1) == 1
        cor = input(i,3:end);
    elseif input(i,1) == 0
        incor = input(i,3:end);
    end

    if i <= round(height(input)/2) && input(i,1) == 1
        correct_firsthalf = [correct_firsthalf;cor];
    elseif i <= round(height(input)/2) && input(i,1) == 0
        incorrect_firsthalf = [incorrect_firsthalf;incor];
    elseif i > round(height(input)/2) && input(i,1) == 1
        correct_secondhalf = [correct_secondhalf;cor];
    elseif i > round(height(input)/2) && input(i,1) == 0
        incorrect_secondhalf = [incorrect_secondhalf;incor];
    end
end

mean_cor_first = mean(correct_firsthalf,1);
mean_cor_second = mean(correct_secondhalf,1);
mean_inc_first = mean(incorrect_firsthalf,1);
mean_inc_second = mean(incorrect_secondhalf,1);
        