function [stat1,stat2] = NERD_STATS(toc,numberoffiles);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
performance = toc;
avg_per_file = performance / numberoffiles;
disp("NERD STATS")
stat1 = fprintf("run time: %f s\n",performance);
stat2 = fprintf("runtime average per file: %f s\n",avg_per_file);

end