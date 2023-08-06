clear all; close all;

load("prl_stream_analysis.mat")
time = prl_stream_analysis.info.time; 


auc = [];
for ii = 1:height(prl_stream_analysis.acq2.cue)
    stream = table2array(prl_stream_analysis.acq2.cue(ii,4:end));
    area = abs(trapz(time,stream));
    auc = [auc; area];
end
auc = array2table(auc, 'VariableNames',{'Cue'});
prl_stream_analysis.acq2.auc(:,'Cue') = auc;

auc = [];
for ii = 1:height(prl_stream_analysis.rev1.cue)
    stream = table2array(prl_stream_analysis.rev1.cue(ii,4:end));
    area = abs(trapz(time,stream));
    auc = [auc; area];
end
auc = array2table(auc, 'VariableNames',{'Cue'});
prl_stream_analysis.rev1.auc(:,'Cue') = auc;






