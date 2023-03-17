clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
warning

BLOCKPATH = 'Z:\DA_Cbrain\WIN115-FR1-5_WIN108-FR1-5';
EXCELPATH = 'Z:\DA_Cbrain\WIN_Cbrain.xlsx';
%Specifiy the correct sheet number and range of data. Make sure the
%timestamps are vertical in the first column
SHEET = 1;
RANGE = "AW1:AW17";
mouse_sex = 0; %1=female, 0=male
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'streams'});
exceldata = readcell(EXCELPATH, 'Sheet', SHEET, 'Range', RANGE);
infusion_ts = cell2mat(exceldata(3:end-2, 1));
end_ts = height(infusion_ts);
% time bins to look for peaks in seconds %
beforeINF = 12;
afterINF = 30;
%Cbrain Variables%
if mouse_sex == 1
    a = 9.41;% for females
elseif mouse_sex == 0
    a = 9.65;% for males
end
alpha = 0.097;
beta = 0.642;

%stretches the infusion timstamps out in K second intervals%
K = 30;
j = 1;

t = infusion_ts(1,1);
u = infusion_ts(end_ts,1);
ts_new = double(t:K:u).';
for i = 1:end_ts
    if i == end_ts
        break
    end
    t = infusion_ts(i,1);
    u = infusion_ts(i+1,1);
    ts_new2 = double(t:K:u).';
    len = height(ts_new2);
    ts_stretch(j:j+len-1) = ts_new2;
    ts_final = ts_stretch.';
    j = j+len;
end
end_ts2 = height(ts_new);
end_ts3 = height(ts_final);

if isfield(data.streams, 'x465C')
    %Stream Stores%
    ISOS1 = 'x405A'; % name of the 405A store
    DOPE1 = 'x465A'; % name of the 465A store
    ISOS2 = 'x405C'; % name of the 405C store
    DOPE2 = 'x465C'; % name of the 465C store
elseif isfield(data.streams, 'b2_B')
    %Stream Stores%
    ISOS1 = 'x405A'; % name of the 405A store
    DOPE1 = 'x470A'; % name of the 465A store
    ISOS2 = 'v2_B'; % name of the 405C store
    DOPE2 = 'b2_B'; % name of the 465C store
elseif isfield(data.streams, 'B2_B')
    %Stream Stores%
    ISOS1 = 'x405A'; % name of the 405A store
    DOPE1 = 'x470A'; % name of the 465A store
    ISOS2 = 'V2_B'; % name of the 405C store
    DOPE2 = 'B2_B'; % name of the 465C store
end

       
N = 100; %Downsample N times
%time array used for all streams%
time = (1:length(data.streams.(DOPE1).data))/data.streams.(DOPE1).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%
t = 30; % time threshold below which we will discard
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(DOPE1).data = data.streams.(DOPE1).data(ind:end);
data.streams.(ISOS1).data = data.streams.(ISOS1).data(ind:end);
data.streams.(DOPE2).data = data.streams.(DOPE2).data(ind:end);
data.streams.(ISOS2).data = data.streams.(ISOS2).data(ind:end);

%downsample streams and time array by N times%
data.streams.(ISOS1).data = downsample(data.streams.(ISOS1).data, N);
data.streams.(DOPE1).data = downsample(data.streams.(DOPE1).data, N);
data.streams.(ISOS2).data = downsample(data.streams.(ISOS2).data, N);
data.streams.(DOPE2).data = downsample(data.streams.(DOPE2).data, N);
time = downsample(time, N);

% time = time(1, 1:end-2);
% data.streams.(DOPE1).data = data.streams.(DOPE1).data(1:end-2);
% data.streams.(DOPE2).data = data.streams.(DOPE2).data(1:end-2);
% % data.streams.(ISOS1).data = data.streams.(ISOS1).data(1:end-1);
% data.streams.(ISOS2).data = data.streams.(ISOS2).data(1:end-2);

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

z465A = zscore(detrend_465A);
z465C = zscore(detrend_465C);

%calculates and plots median absolute deviation for both 465 signals%
MAD1 = mad(detrend_465A, 1);
MAD2 = mad(detrend_465C, 1);
MAD3 = mad(z465A, 1);
MAD4 = mad(z465C, 1);
MeanMAD = mad(z465A);
[pks,locs,w,p] = findpeaks(z465A, time, "MinPeakHeight", MAD3);
[pks2,locs2,w2,p2] = findpeaks(z465C, time, 'MinPeakHeight', MAD4);

%plots entire streams with peak indicators% 
subplot(2,1,1);
findpeaks(z465A, time, "MinPeakHeight", MAD3);
subplot(2,1,2)
findpeaks(z465C, time, 'MinPeakHeight', MAD4);

%loops through timestamps in excel sheet to find peaks between each
%timestamp segment%
for x = 2:end_ts
    if x == end_ts
        break
    end
    x1 = infusion_ts(x,1);
    x2 = infusion_ts(x+1,1);
    %only use the if statement below for PR animals%
    if x2 - x1 < 30
        continue
    end
    x3 = infusion_ts(x,1)-beforeINF;
    x4 = infusion_ts(x,1);
    x5 = infusion_ts(x,1);
    x6 = infusion_ts(x,1)+afterINF;
    [c, index1] = min(abs(time-x1));
    [c2, index2] = min(abs(time-x2));
    [c3, index3] = min(abs(time-x3));
    [c4, index4] = min(abs(time-x4));
    [c5, index5] = min(abs(time-x5));
    [c6, index6] = min(abs(time-x6));
    DOPE1_sig = z465A(1,index1:index2);
    DOPE2_sig = z465C(1,index1:index2);
    DOPE1_time = time(1,index1:index2);
    DOPE3_sig = z465A(1,index3:index4);
    DOPE4_sig = z465C(1,index3:index4);
    DOPE2_time = time(1,index3:index4);
    DOPE5_sig = z465A(1,index5:index6);
    DOPE6_sig = z465C(1,index5:index6);
    DOPE3_time = time(1,index5:index6);
%     calculates and plots median absolute deviation for both 465 signals%
%     [pks3, locs3] = findpeaks(DOPE1_sig, DOPE1_time, 'MinPeakHeight', MAD3);
%     [pks4, locs4] = findpeaks(DOPE2_sig, DOPE1_time, 'MinPeakHeight', MAD4);
%     [pks5, locs5] = findpeaks(DOPE3_sig, DOPE2_time, 'MinPeakHeight', MAD3);
%     [pks6, locs6] = findpeaks(DOPE4_sig, DOPE2_time, 'MinPeakHeight', MAD4);
%     [pks7, locs7] = findpeaks(DOPE5_sig, DOPE3_time, 'MinPeakHeight', MAD3);
%     [pks8, locs8] = findpeaks(DOPE6_sig, DOPE3_time, 'MinPeakHeight', MAD4);
%     
%     DOPE1_numpeak = length(pks3);
%     DOPE1_pks_min = (DOPE1_numpeak/(x2-x1))*60;
%     DOPE2_numpeak = length(pks4);
%     DOPE2_pks_min = (DOPE2_numpeak/(x2-x1))*60;
%     DOPE3_numpeak = length(pks5);
%     DOPE3_pks_min = (DOPE3_numpeak/(x4-x3))*60;
%     DOPE4_numpeak = length(pks6);
%     DOPE4_pks_min = (DOPE4_numpeak/(x4-x3))*60;
%     DOPE5_numpeak = length(pks7);
%     DOPE5_pks_min = (DOPE5_numpeak/(x6-x5))*60;
%     DOPE6_numpeak = length(pks8);
%     DOPE6_pks_min = (DOPE6_numpeak/(x6-x5))*60;
%     
%     peak_analysis(x,:) = [x1 DOPE1_pks_min x3 DOPE3_pks_min x6,DOPE5_pks_min ...
%         x1 DOPE2_pks_min x3 DOPE4_pks_min x6 DOPE6_pks_min];
    
end
Cbrain = infusion_ts/60;
length_Cbrain = length(Cbrain);


for j = 2:length_Cbrain
    if j == length_Cbrain
        break
    end
    m = Cbrain(j,:);
    m1 = Cbrain(j+1,:);
    Csum(j) = (m1-m); 
    
end
Csum = Csum.';
for y = 1:height(Csum)
    diff1 = Csum(y,1);
    diff2 = Csum(y,1) - (12/60);
    diff3 = Csum(y,1) + (30/60);
    XX = a*(exp((-alpha)*diff1)-exp((-beta)*diff1));
    YY = a*(exp((-alpha)*diff2)-exp((-beta)*diff2));
    ZZ = a*(exp((-alpha)*diff3)-exp((-beta)*diff3));
    if y > 1
        XX = XX + Csum(y-1,2);
        YY = YY + Csum(y-1,3);
        ZZ = ZZ + Csum(y-1,4);
    end
    Csum(y,2) = XX;
    Csum(y,3) = YY;
    Csum(y,4) = ZZ;
end
Csum(1,3) = 0;

% Plots DA in 30 second increments %
for x = 1:end_ts3
    if x == end_ts3-1
        break
    end
    x1 = ts_final(x,1);
    [c, index] = min(abs(time-x1));
    
    DOPE_z_loc = z465A(1,index);
    DOPE_z_loc2 = z465C(1,index);
    DOPE_dFF_loc = detrend_465A(1,index);
    DOPE_dFF_loc2 = detrend_465C(1,index);
    DOPE_time_m = ts_final(x,1)/60;
    DOPE_time_s = ts_final(x,1);
    dope_stream(x,:) = [DOPE_time_s DOPE_time_m DOPE_z_loc DOPE_dFF_loc DOPE_z_loc2 DOPE_dFF_loc2];

end

