clear all; clc; close all;

BLOCKPATH = '/Users/brandon/Desktop/DA PRL/Tanks/Condensed PRL (pilot)/DA19_Session2_DA20_Session2_JZL';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
box_number = 3;
training_logic = 2; %Training = 1, all else = 2%
st1 = cellstr(sprintf('%s','St1/'));
crewa = cellstr(sprintf('%s','cRewA'));
cnorewa = cellstr(sprintf('%s','cNoRewA'));
st2 = cellstr(sprintf('%s','St2/'));
crewc = cellstr(sprintf('%s','cRewC'));
cnorewc = cellstr(sprintf('%s','cNoRewC'));
pe1 = cellstr(sprintf('%s','Pe1/'));
pe2 = cellstr(sprintf('%s','Pe2/'));
cl1 = cellstr(sprintf('%s','CL1/'));
cl2 = cellstr(sprintf('%s','CL2/'));

epoch_listA = [st1 crewa cnorewa];
epoch_listC = [st2 crewc cnorewc];
training_epocA = [st1 pe1 cl1];
training_epocC = [st2 pe2 cl2];

%Makes separate epocs for different trial outcomes
correct_ts_A = data.epocs.CL1_.onset;
pellet_ts_A = data.epocs.Pe1_.onset;
correct_rewardedA = zeros(double(height(pellet_ts_A)));
correct_norewardA = zeros(double(height(pellet_ts_A)));
for i = 1:height(correct_ts_A)
    var1 = correct_ts_A(i,:);
    var2 = pellet_ts_A;
    [f, ia] = ismember(var1, pellet_ts_A);
    if f == 1
        correct_rewardedA(i,:) = var1;
    elseif f < 1
        correct_norewardA(i,:) = var1;
    end  
end
correct_rewardedA = nonzeros(correct_rewardedA(:,1));
correct_norewardA = nonzeros(correct_norewardA(:,1));
data.epocs.cRewA.name = 'cRewA';
data.epocs.cRewA.onset = correct_rewardedA;
data.epocs.cRewA.offset = correct_rewardedA+1;
data.epocs.cRewA.data = ones(height(correct_rewardedA));
data.epocs.cNoRewA.name = 'cNoRewA';
data.epocs.cNoRewA.onset = correct_norewardA;
data.epocs.cNoRewA.offset = correct_norewardA+1;
data.epocs.cNoRewA.data = ones(height(correct_norewardA))*2;

correct_ts_C = data.epocs.CL2_.onset;
pellet_ts_C = data.epocs.Pe2_.onset;
correct_rewardedC = zeros(double(height(pellet_ts_C)));
correct_norewardC = zeros(double(height(pellet_ts_C)));
for i = 1:height(correct_ts_C)
    var3 = correct_ts_C(i,:);
    var4 = pellet_ts_C;
    [g, ia] = ismember(var3, pellet_ts_C);
    if g == 1
        correct_rewardedC(i,:) = var3;
    elseif g < 1
        correct_norewardC(i,:) = var3;
    end  
end
correct_rewardedC = nonzeros(correct_rewardedC(:,1));
correct_norewardC = nonzeros(correct_norewardC(:,1));
data.epocs.cRewC.name = 'cRewC';
data.epocs.cRewC.onset = correct_rewardedC;
data.epocs.cRewC.offset = correct_rewardedC+1;
data.epocs.cRewC.data = ones(height(correct_rewardedC));
data.epocs.cNoRewC.name = 'cNoRewC';
data.epocs.cNoRewC.onset = correct_norewardC;
data.epocs.cNoRewC.offset = correct_norewardC+1;
data.epocs.cNoRewC.data = ones(height(correct_norewardC))*2;

for d = 1:3
    if box_number == 3 && training_logic == 2
        STREAM_STORE1 = 'x405A'; % name of the 405 store
        STREAM_STORE2 = 'x465A'; % name of the 465 store
        d1 = epoch_listA(:,d);
    elseif box_number == 4 && training_logic == 2
        STREAM_STORE1 = 'x405C'; % name of the 405 store
        STREAM_STORE2 = 'x465C'; % name of the 465 store
        d1 = epoch_listC(:,d);
    elseif box_number == 3 && training_logic == 1
        STREAM_STORE1 = 'x405A'; % name of the 405 store
        STREAM_STORE2 = 'x465A'; % name of the 465 store
        d1 = training_epocA(:,d);
    elseif box_number == 4 && training_logic == 1
        STREAM_STORE1 = 'x405C'; % name of the 405 store
        STREAM_STORE2 = 'x465C'; % name of the 465 store
        d1 = training_epocC(:,d);
    end
    REF_EPOC = char(d1); % Stimulation event to center on
    TRANGE = [-2 9]; %window size [start time relative to epoc onset, entire duration]
    ARANGE = [1 1];
    BASELINE_PER = [-5 -1]; % baseline period before stim
    ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
    ARTIFACT465 = Inf;% variable created for artifact removal for 465 store
    % Use TDTfilter to extract data around our epoc event
    % Using the 'TIME' parameter extracts data only from the time range around
    % our epoc event. Use the 'VALUES' parameter to specify allowed values of
    % the REF_EPOC to extract.  For stream events, the chunks of data are 
    % stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
    data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);  
   
    
    % Optionally remove artifacts. If any waveform is above ARTIFACT level, or
    % below -ARTIFACT level, remove it from the data set.
    art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
    art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
    good = ~art1 & ~art2;
    data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);
    
    art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT465), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
    art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT465), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
    good2 = ~art1 & ~art2;
    data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);
    
    numArtifacts = sum(~good) + sum(~good2);
    
    %%
    % Applying a time filter to a uniformly sampled signal means that the
    % length of each segment could vary by one sample.  Let's find the minimum
    % length so we can trim the excess off before calculating the mean.
    minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
    minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
    data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
    data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);
    
    allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');
    
    % downsample 10x and average 405 signal
    N = 10;
    F405 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength1 = size(F405,2);
    
    % Create mean signal, standard error of signal, and DC offset of 405 signal
    meanSignal1 = mean(F405);
    stdSignal1 = std(double(F405))/sqrt(size(F405,1));
    dcSignal1 = mean(meanSignal1);
    
    % downsample 10x and average 465 signal
    allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
    F465 = zeros(size(allSignals(:,1:N:end-N+1)));
    for ii = 1:size(allSignals,1)
        F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
    end
    minLength2 = size(F465,2);
    
    % Create mean signal, standard error of signal, and DC offset of 465 signal
    meanSignal2 = mean(F465);
    stdSignal2 = std(double(F465))/sqrt(size(F465,1));
    dcSignal2 = mean(meanSignal2);
    
    %% Plot Epoch Averaged Response
    
    % Create the time vector for each stream store
    ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
    ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;
    
    % Subtract DC offset to get signals on top of one another
    meanSignal1 = meanSignal1 - dcSignal1;
    meanSignal2 = meanSignal2 - dcSignal2;
    
    bls = polyfit(F465(1:end), F405(1:end), 1);
    Y_fit_all = bls(1) .* F405 + bls(2);
    Y_dF_all = F465 - Y_fit_all;
    
    zall = zeros(size(Y_dF_all));
    for i = 1:size(Y_dF_all,1)
        ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
        zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
        zsd = std(Y_dF_all(i,ind)); % baseline period stdev
        zall(i,:)=(Y_dF_all(i,:) - zb)/zsd; % Z score per bin
    end
    
    % Standard error of the z-score
    meanZall = mean(zall);
    zerror = std(zall)/sqrt(size(zall,1));
    
    if box_number == 3
        mean_epoch = meanSignal1;
    elseif box_number == 4
        mean_epoch = meanSignal2;
    end
    mean_epoch_analysis(:,d) = mean_epoch;
end
mean_epoch_analysis = array2table(mean_epoch_analysis, 'VariableNames', {'Cue', 'CorRew', 'CorNoRew'});
