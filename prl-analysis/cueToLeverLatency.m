correct = data.epocs.CL1_.onset;
incorrect = data.epocs.IL1_.onset;
cue = data.epocs.St1_.onset;
levers = [correct;incorrect];
levers = sort(levers);
latency = [];
for i = 5:height(cue)
    x = levers(i,1) - cue(i,1);
    latency = [latency;x];
end
    