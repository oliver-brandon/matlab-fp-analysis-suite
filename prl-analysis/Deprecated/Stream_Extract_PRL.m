clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
warning
%Windows%
BLOCKPATH = '/Users/brandon/Desktop/DA_PRL/Tanks/DA67_Rev1_JZL_DA69_Rev1_JZL';
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box_number = 4; % 3 = FibPho1, 4 = FibPho2
training_logic = 2; % Training = 1, all else (acq/rev) = 2%
timeWindow = 5; % the number of seconds after the onset of a TTL to analyze
baseline = 2; % baseline signal to include before TTL 
baselineZ = [5 1];
N = 100; %Downsample N times
minArrayLen = 72; 
%array column length definition to eliminate error produced
%when trying to fill array with stream snips of different lengths 
%(negative relationship with N (downsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});


ISOS1 = 'x405A';
DOPE1 = 'x465A';
ISOS2 = 'x405C';
DOPE2 = 'x465C';
% data.streams.(DOPE2).data = data.streams.(DOPE2).data(1:length(data.streams.(ISOS2).data));
master_amp_analysisA = zeros(30,5);
master_amp_analysisC = zeros(30,5);
cRew1AMP = [];
cNoRew1AMP = [];
iRew1AMP = [];
iNoRew1AMP = [];
cRew2AMP = [];
cNoRew2AMP = [];
iRew2AMP = [];
iNoRew2AMP = [];
cRew1STREAM = [];
cNoRew1STREAM = [];
iRew1STREAM = [];
iNoRew1STREAM = [];
cRew2STREAM = [];
cNoRew2STREAM = [];
iRew2STREAM = [];
iNoRew2STREAM = [];
if box_number == 3 && training_logic == 2
    cueTSA = data.epocs.St1_.onset;
    %Makes separate epocs for different trial outcomes
    if isfield(data.epocs, 'CL1_') == 0
        correct_ts_A = 0;
    else
        correct_ts_A = data.epocs.CL1_.onset;
    end
    
    if isfield(data.epocs, 'IL1_') == 0
        incorrect_ts_A = 0;
    else
        incorrect_ts_A = data.epocs.IL1_.onset;
    end
    pellet_ts_A = data.epocs.Pe1_.onset;
    correct_rewardedA = zeros(double(height(pellet_ts_A)));
    correct_norewardA = zeros(double(height(pellet_ts_A)));
    incorrect_rewardedA = zeros(double(height(pellet_ts_A)));
    incorrect_norewardA = zeros(double(height(pellet_ts_A)));
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
    for i = 1:height(incorrect_ts_A)
        var1 = incorrect_ts_A(i,:);
        var2 = pellet_ts_A;
        [f, ia] = ismember(var1, pellet_ts_A);
        if f == 1
            incorrect_rewardedA(i,:) = var1;
        elseif f < 1
            incorrect_norewardA(i,:) = var1;
        end
    end
    
    correct_rewardedA = nonzeros(correct_rewardedA(:,1));
    correct_norewardA = nonzeros(correct_norewardA(:,1));
    incorrect_rewardedA = nonzeros(incorrect_rewardedA(:,1));
    incorrect_norewardA = nonzeros(incorrect_norewardA(:,1));

    if isempty(correct_rewardedA) == 1
    correct_rewardedA = 0;
    end
    if isempty(correct_norewardA) == 1
        correct_norewardA = 0;
    end
    if isempty(incorrect_rewardedA) == 1
        incorrect_rewardedA = 0;
    end
    if isempty(incorrect_norewardA) == 1
        incorrect_norewardA = 0;
    end
    data.epocs.cRewA.name = 'cRewA';
    data.epocs.cRewA.onset = correct_rewardedA;
    data.epocs.cRewA.offset = correct_rewardedA+1;
    data.epocs.cRewA.data = ones(height(correct_rewardedA));
    data.epocs.cNoRewA.name = 'cNoRewA';
    data.epocs.cNoRewA.onset = correct_norewardA;
    data.epocs.cNoRewA.offset = correct_norewardA+1;
    data.epocs.cNoRewA.data = ones(height(correct_norewardA))*2;
    data.epocs.iRewA.name = 'iRewA';
    data.epocs.iRewA.onset = incorrect_rewardedA;
    data.epocs.iRewA.offset = incorrect_rewardedA+1;
    data.epocs.iRewA.data = ones(height(incorrect_rewardedA))*3;
    data.epocs.iNoRewA.name = 'iNoRewA';
    data.epocs.iNoRewA.onset = incorrect_norewardA;
    data.epocs.iNoRewA.offset = incorrect_norewardA+1;
    data.epocs.iNoRewA.data = ones(height(incorrect_norewardA))*4;
elseif box_number == 4 && training_logic == 2
    cueTSC = data.epocs.St2_.onset;
    %Makes separate epocs for different trial outcomes
    if isfield(data.epocs, 'CL2_') == 0
        correct_ts_C = 0;
    else
        correct_ts_C = data.epocs.CL2_.onset;
    end
    
    if isfield(data.epocs, 'IL2_') == 0
        incorrect_ts_C = 0;
    else
        incorrect_ts_C = data.epocs.IL2_.onset;
    end
    pellet_ts_C = data.epocs.Pe2_.onset;
    correct_rewardedC = zeros(double(height(pellet_ts_C)));
    correct_norewardC = zeros(double(height(pellet_ts_C)));
    incorrect_rewardedC = zeros(double(height(pellet_ts_C)));
    incorrect_norewardC = zeros(double(height(pellet_ts_C)));
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
    for i = 1:height(incorrect_ts_C)
        var3 = incorrect_ts_C(i,:);
        var4 = pellet_ts_C;
        [g, ia] = ismember(var3, pellet_ts_C);
        if g == 1
            incorrect_rewardedC(i,:) = var3;
        elseif g < 1
            incorrect_norewardC(i,:) = var3;
        end
    end
    correct_rewardedC = nonzeros(correct_rewardedC(:,1));
    correct_norewardC = nonzeros(correct_norewardC(:,1));
    incorrect_rewardedC = nonzeros(incorrect_rewardedC(:,1));
    incorrect_norewardC = nonzeros(incorrect_norewardC(:,1));
    if isempty(correct_rewardedC) == 1
        correct_rewardedC = 0;
    end
    if isempty(correct_norewardC) == 1
        correct_norewardC = 0;
    end
    if isempty(incorrect_rewardedC) == 1
        incorrect_rewardedC = 0;
    end
    if isempty(incorrect_norewardC) == 1
        incorrect_norewardC = 0;
    end
    data.epocs.cRewC.name = 'cRewC';
    data.epocs.cRewC.onset = correct_rewardedC;
    data.epocs.cRewC.offset = correct_rewardedC+1;
    data.epocs.cRewC.data = ones(height(correct_rewardedC));
    data.epocs.cNoRewC.name = 'cNoRewC';
    data.epocs.cNoRewC.onset = correct_norewardC;
    data.epocs.cNoRewC.offset = correct_norewardC+1;
    data.epocs.cNoRewC.data = ones(height(correct_norewardC))*2;
    data.epocs.iRewC.name = 'iRewC';
    data.epocs.iRewC.onset = incorrect_rewardedC;
    data.epocs.iRewC.offset = incorrect_rewardedC+1;
    data.epocs.iRewC.data = ones(height(incorrect_rewardedC))*3;
    data.epocs.iNoRewC.name = 'iNoRewC';
    data.epocs.iNoRewC.onset = incorrect_norewardC;
    data.epocs.iNoRewC.offset = incorrect_norewardC+1;
    data.epocs.iNoRewC.data = ones(height(incorrect_norewardC))*4;
end



%time array used for all streams%
time1 = (1:length(data.streams.(DOPE1).data))/data.streams.(DOPE1).fs;
time2 = (1:length(data.streams.(DOPE2).data))/data.streams.(DOPE2).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%
t = 0; % time threshold below which we will discard
ind = find(time1>t,1);% find first index of when time crosses threshold
ind2 = find(time2>t,1);
time1 = time1(ind:end); % reformat vector to only include allowed time
time2 = time2(ind2:end);
data.streams.(DOPE1).data = data.streams.(DOPE1).data(ind:end);
data.streams.(ISOS1).data = data.streams.(ISOS1).data(ind:end);
data.streams.(DOPE2).data = data.streams.(DOPE2).data(ind2:end);
data.streams.(ISOS2).data = data.streams.(ISOS2).data(ind2:end);

%downsample streams and time array by N times%
data.streams.(ISOS1).data = downsample(data.streams.(ISOS1).data, N);
data.streams.(DOPE1).data = downsample(data.streams.(DOPE1).data, N);
data.streams.(ISOS2).data = downsample(data.streams.(ISOS2).data, N);
data.streams.(DOPE2).data = downsample(data.streams.(DOPE2).data, N);
time1 = downsample(time1, N);
time2 = downsample(time2, N);
% start_val = -2;
% inc = (1:minArrayLen) / data.streams.(DOPE1).fs*N;
% stop_val = (minArrayLen-1)*inc + start_val;
% v = [start_val:inc:stop_val];
ts1 = -baseline + (1:minArrayLen) / data.streams.(DOPE1).fs*N;
ts2 = -baseline + (1:minArrayLen) / data.streams.(DOPE2).fs*N;
%detrend & dFF%
%465A%
bls = polyfit(data.streams.(ISOS1).data,data.streams.(DOPE1).data,1);
Y_fit_all = bls(1) .* data.streams.(ISOS1).data + bls(2);
Y_dF_all = data.streams.(DOPE1).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
detrend_465A = detrend(dFF);

%465C%
bls2 = polyfit(data.streams.(ISOS2).data,data.streams.(DOPE2).data,1);
Y_fit_all2 = bls2(1) .* data.streams.(ISOS2).data + bls2(2);
Y_dF_all2 = data.streams.(DOPE2).data - Y_fit_all2; %dF (units mV) is not dFF
dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
std_dFF2 = std(double(dFF2));
detrend_465C = detrend(dFF2);
%Makes streams/time same length
% minStreamLen = min(numel(detrend_465A), numel(detrend_465C));
% time1 = time1(1:minStreamLen);
% time2 = time1;
% detrend_465A = detrend_465A(1:minStreamLen);
% detrend_465C = detrend_465C(1:minStreamLen);

z465A = zscore(data.streams.(DOPE1).data);
z465C = zscore(data.streams.(DOPE2).data);

%calculates and plots median absolute deviation for both 465 signals%
MAD1 = mad(z465A, 1);
MAD2 = mad(z465C, 1);

for d = 1:5
    if box_number == 3 && training_logic == 2
        cue1STREAM = [];
        for e = 1:height(cueTSA)
            e1 = data.epocs.St1_.onset(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = data.epocs.St1_.onset(e,1)-baselineZ(:,2);
            e4 = data.epocs.St1_.onset(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            DOPE1_zBase = detrend_465A(1,ind4:ind3);
            DOPE1_signal = detrend_465A(1,ind1:ind2);
            DOPE1_time = time1(1,ind1:ind2);
            zb = mean(DOPE1_zBase);
            zsd = std(DOPE1_zBase);
            zfinal = (DOPE1_signal - zb)/zsd;
            cue1STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE1_signal, DOPE1_time, 'MinPeakHeight', MAD1);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE1_signal = zscore(DOPE1_signal);
            cue1AMP(e,:) = max(zfinal);
        end
        for e = 1:height(correct_rewardedA)
            e1 = correct_rewardedA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_rewardedA(e,1)-baselineZ(:,2);
            e4 = correct_rewardedA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            DOPE1_zBase = detrend_465A(1,ind4:ind3);
            DOPE1_signal = detrend_465A(1,ind1:ind2);
            DOPE1_time = time1(1,ind1:ind2);
            zb = mean(DOPE1_zBase);
            zsd = std(DOPE1_zBase);
            zfinal = (DOPE1_signal - zb)/zsd;
            if correct_rewardedA == 0
                break
            end
            cRew1STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE1_signal, DOPE1_time, 'MinPeakHeight', MAD1);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE1_signal = zscore(DOPE1_signal);
            cRew1AMP(e,:) = max(zfinal);
        end
        for e = 1:height(correct_norewardA)
            e1 = correct_norewardA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_norewardA(e,1)-baselineZ(:,2);
            e4 = correct_norewardA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            DOPE1_zBase = detrend_465A(1,ind4:ind3);
            DOPE1_signal = detrend_465A(1,ind1:ind2);
            DOPE1_time = time1(1,ind1:ind2);
            zb = mean(DOPE1_zBase);
            zsd = std(DOPE1_zBase);
            zfinal = (DOPE1_signal - zb)/zsd;
            if correct_norewardA == 0
                break
            end
            cNoRew1STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE1_signal, DOPE1_time, 'MinPeakHeight', MAD1);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE1_signal = zscore(DOPE1_signal);
            cNoRew1AMP(e,:) = max(zfinal);
        end
        for e = 1:height(incorrect_rewardedA)
            e1 = incorrect_rewardedA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_rewardedA(e,1)-baselineZ(:,2);
            e4 = incorrect_rewardedA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            DOPE1_zBase = detrend_465A(1,ind4:ind3);
            DOPE1_signal = detrend_465A(1,ind1:ind2);
            DOPE1_time = time1(1,ind1:ind2);
            zb = mean(DOPE1_zBase);
            zsd = std(DOPE1_zBase);
            zfinal = (DOPE1_signal - zb)/zsd;
            if incorrect_rewardedA == 0
                break
            end
            iRew1STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE1_signal, DOPE1_time, 'MinPeakHeight', MAD1);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE1_signal = zscore(DOPE1_signal);
            iRew1AMP(e,:) = max(zfinal);
        end
        for e = 1:height(incorrect_norewardA)
            e1 = incorrect_norewardA(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_norewardA(e,1)-baselineZ(:,2);
            e4 = incorrect_norewardA(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time1-e1));
            [c2,ind2] = min(abs(time1-e2));
            [c3,ind3] = min(abs(time1-e3));
            [c4,ind4] = min(abs(time1-e4));
            DOPE1_zBase = detrend_465A(1,ind4:ind3);
            DOPE1_signal = detrend_465A(1,ind1:ind2);
            DOPE1_time = time1(1,ind1:ind2);
            zb = mean(DOPE1_zBase);
            zsd = std(DOPE1_zBase);
            zfinal = (DOPE1_signal - zb)/zsd;
            iNoRew1STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE1_signal, DOPE1_time, 'MinPeakHeight', MAD1);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE1_signal = zscore(DOPE1_signal);
            iNoRew1AMP(e,:) = max(zfinal);
        end
    elseif box_number == 4 && training_logic == 2
        for e = 1:height(cueTSC)
            e1 = data.epocs.St2_.onset(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            e3 = data.epocs.St2_.onset(e,1)-baselineZ(:,2);
            e4 = data.epocs.St2_.onset(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time2-e1));
            [c2,ind2] = min(abs(time2-e2));
            [c3,ind3] = min(abs(time2-e3));
            [c4,ind4] = min(abs(time2-e4));
            DOPE2_zBase = detrend_465C(1,ind4:ind3);
            DOPE2_signal = detrend_465C(1,ind1:ind2);
            DOPE2_time = time2(1,ind1:ind2);
            zb = mean(DOPE2_zBase);
            zsd = std(DOPE2_zBase);
            zfinal = (DOPE2_signal - zb)/zsd;
            cue2STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE2_signal, DOPE2_time, 'MinPeakHeight', MAD2);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE2_signal = zscore(DOPE2_signal);
            cue2AMP(e,:) = max(zfinal);
        end
        for e = 1:height(correct_rewardedC)
            e1 = correct_rewardedC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_rewardedC(e,1)-baselineZ(:,2);
            e4 = correct_rewardedC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time2-e1));
            [c2,ind2] = min(abs(time2-e2));
            [c3,ind3] = min(abs(time2-e3));
            [c4,ind4] = min(abs(time2-e4));
            DOPE2_zBase = detrend_465C(1,ind4:ind3);
            DOPE2_signal = detrend_465C(1,ind1:ind2);
            DOPE2_time = time2(1,ind1:ind2);
            zb = mean(DOPE2_zBase);
            zsd = std(DOPE2_zBase);
            zfinal = (DOPE2_signal - zb)/zsd;
            if correct_rewardedC == 0
                break
            end
            cRew2STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE2_signal, DOPE2_time, 'MinPeakHeight', MAD2);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE2_signal = zscore(DOPE2_signal);
            cRew2AMP(e,:) = max(zfinal);
        end
        for e = 1:height(correct_norewardC)
            e1 = correct_norewardC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = correct_norewardC(e,1)-baselineZ(:,2);
            e4 = correct_norewardC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time2-e1));
            [c2,ind2] = min(abs(time2-e2));
            [c3,ind3] = min(abs(time2-e3));
            [c4,ind4] = min(abs(time2-e4));
            DOPE2_zBase = detrend_465C(1,ind4:ind3);
            DOPE2_signal = detrend_465C(1,ind1:ind2);
            DOPE2_time = time2(1,ind1:ind2);
            zb = mean(DOPE2_zBase);
            zsd = std(DOPE2_zBase);
            zfinal = (DOPE2_signal - zb)/zsd;
            if correct_norewardC == 0
                break
            end
            cNoRew2STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE2_signal, DOPE2_time, 'MinPeakHeight', MAD2);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE2_signal = zscore(DOPE2_signal);
            cNoRew2AMP(e,:) = max(zfinal);
        end
        for e = 1:height(incorrect_rewardedC)
            e1 = incorrect_rewardedC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_rewardedC(e,1)-baselineZ(:,2);
            e4 = incorrect_rewardedC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time2-e1));
            [c2,ind2] = min(abs(time2-e2));
            [c3,ind3] = min(abs(time2-e3));
            [c4,ind4] = min(abs(time2-e4));
            DOPE2_zBase = detrend_465C(1,ind4:ind3);
            DOPE2_signal = detrend_465C(1,ind1:ind2);
            DOPE2_time = time2(1,ind1:ind2);
            zb = mean(DOPE2_zBase);
            zsd = std(DOPE2_zBase);
            zfinal = (DOPE2_signal - zb)/zsd;
            if incorrect_rewardedC == 0
                break
            end
            iRew2STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE2_signal, DOPE2_time, 'MinPeakHeight', MAD2);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE2_signal = zscore(DOPE2_signal);
            iRew2AMP(e,:) = max(zfinal);
        end
        for e = 1:height(incorrect_norewardC)
            e1 = incorrect_norewardC(e,1)-baseline;
            e2 = e1+timeWindow+baseline;
            if e1 == 0
                continue
            end
            e3 = incorrect_norewardC(e,1)-baselineZ(:,2);
            e4 = incorrect_norewardC(e,1)-baselineZ(:,1);
            [c1,ind1] = min(abs(time2-e1));
            [c2,ind2] = min(abs(time2-e2));
            [c3,ind3] = min(abs(time2-e3));
            [c4,ind4] = min(abs(time2-e4));
            DOPE2_zBase = detrend_465C(1,ind4:ind3);
            DOPE2_signal = detrend_465C(1,ind1:ind2);
            DOPE2_time = time2(1,ind1:ind2);
            zb = mean(DOPE2_zBase);
            zsd = std(DOPE2_zBase);
            zfinal = (DOPE2_signal - zb)/zsd;
            if incorrect_norewardC == 0
                break
            end
            iNoRew2STREAM(e,:) = zfinal(1:minArrayLen);
%             [pks1, locs1] = findpeaks(DOPE2_signal, DOPE2_time, 'MinPeakHeight', MAD2);
%             if isempty(pks1) == 1
%                 pks1 = 0;
%             end
%             zDOPE2_signal = zscore(DOPE2_signal);
            iNoRew2AMP(e,:) = max(zfinal);
        end
%     elseif box_number == 3 && training_logic == 1
%     elseif box_number == 4 && training_logic == 1
    end
    if box_number == 3 && training_logic == 2
        if height(cue1STREAM) > 1
            meanCue1 = mean(cue1STREAM);
        else 
            meanCue1 = cue1STREAM;
        end 
        if height(cRew1STREAM) > 1
            meanCRew1 = mean(cRew1STREAM);
        elseif correct_rewardedA == 0
            meanCRew1 = zeros([1 minArrayLen]);
        else
            meanCRew1 = cRew1STREAM;
        end
        if height(cNoRew1STREAM) > 1
            meanCNRew1 = mean(cNoRew1STREAM);
        elseif correct_norewardA == 0
            meanCNRew1 = zeros([1 minArrayLen]);
        else
            meanCNRew1 = cNoRew1STREAM;
        end
        if height(iRew1STREAM) > 1
            meanIRew1 = mean(iRew1STREAM);
        elseif incorrect_rewardedA == 0
            meanIRew1 = zeros([1 minArrayLen]);
        else 
            meanIRew1 = iRew1STREAM;
        end
        if height(iNoRew1STREAM) > 1
            meanINRew1 = mean(iNoRew1STREAM);
        elseif incorrect_norewardA == 0
            meanINRew1 = zeros([1 minArrayLen]);
        else
            meanINRew1 = iNoRew1STREAM;
        end
        master_amp_analysisA = table(mean(cue1AMP), mean(cRew1AMP), mean(cNoRew1AMP), mean(iRew1AMP), mean(iNoRew1AMP),...
            'VariableNames', {'Cue Amp','cRew Amp','cNoRew Amp','iRew','iNoRew'});
        master_stream_analysisA = table(meanCue1', meanCRew1', meanCNRew1', meanIRew1',...
            meanINRew1', 'VariableNames', {'cue', 'cRew', 'cNoRew', 'iRew', 'iNoRew'});
    elseif box_number == 4 && training_logic == 2
        if height(cue2STREAM) > 1
            meanCue2 = mean(cue2STREAM);
        else 
            meanCue2 = cue2STREAM;
        end 
        if height(cRew2STREAM) > 1
            meanCRew2 = mean(cRew2STREAM);
        elseif correct_rewardedC == 0
            meanCRew2 = zeros([1 minArrayLen]);
        else
            meanCRew2 = cRew2STREAM;
        end
        if height(cNoRew2STREAM) > 1
            meanCNRew2 = mean(cNoRew2STREAM);
        elseif correct_norewardC == 0
            meanCNRew2 = zeros([1 minArrayLen]);
        else
            meanCNRew2 = cNoRew2STREAM;
        end
        if height(iRew2STREAM) > 1
            meanIRew2 = mean(iRew2STREAM);
        elseif incorrect_rewardedC == 0
            meanIRew2 = zeros([1 minArrayLen]);
        else 
            meanIRew2 = iRew2STREAM;
        end
        if height(iNoRew2STREAM) > 1
            meanINRew2 = mean(iNoRew2STREAM);
        elseif incorrect_norewardC == 0
            meanINRew2 = zeros([1 minArrayLen]);
        else
            meanINRew2 = iNoRew2STREAM;
        end
        master_amp_analysisC = table(mean(cue2AMP), mean(cRew2AMP), mean(cNoRew2AMP), mean(iRew2AMP), mean(iNoRew2AMP),...
            'VariableNames', {'Cue Amp','cRew Amp','cNoRew Amp','iRew','iNoRew'});
        master_stream_analysisC = table(meanCue2', meanCRew2', meanCNRew2', meanIRew2',...
            meanINRew2', 'VariableNames', {'cue', 'cRew', 'cNoRew', 'iRew', 'iNoRew'});
    end
end
disp("DONE")
