% Grab_signalTest.m created by Brandon L. Oliver, M.A.
% Used to quickly test Grab sensor flourescence from TDT fiberphotometry tanks.
% Plots the raw, dFF, and detrended dFF on seperate plots and saves a
% figure to a user designated directory. Can also downsample if necessary
% by changing the value of 'N' below.
clear
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figsavetype = '.pdf'; % can change to '.jpg', '.fig', etc.
Grab_Name = 'GrabDA4.4'; % example: 'GrabDA4.4'
t = 5; % first t seconds are discarded to remove laser on artifact
N = 1; % downsample signal N times
ISOS = 'x405C'; % set name of isosbestic signal
Grab = 'x465C'; % set name of Grab signal
x = 60; % start of window (s)
y = 70; % end of window (s)
fontSize = 10; % font size for figure ylabels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gets tank from UI pop-up window
TANK_NAME = uigetdir(pwd, 'Select a tank to plot');
figsavepath = strcat(TANK_NAME,'/');
[~,name,~] = fileparts(TANK_NAME);
TITLE = strrep(name,'_',' ');
data = TDTbin2mat(TANK_NAME, 'TYPE', {'streams'});
ISOS_raw = data.streams.(ISOS).data;
Grab_raw = data.streams.(Grab).data;

% Checks for unequal isosbestic and Grab sensor signal length correcting if
% necessary
if length(Grab_raw) < length(ISOS_raw)
    disp('Isosbestic signal array is longer than Grab signal array')
    ISOS_raw = ISOS_raw(1:length(Grab_raw));
    disp('Corrected.')
elseif length(Grab_raw) > length(ISOS_raw)
    disp('Isosbestic signal array is shorter than Grab signal array')
    Grab_raw = Grab_raw(1:length(ISOS_raw));
    disp('Corrected.')
else
    disp('Isosbestic and Grab signal arrays are of equal size')
    disp('No correction necessary.')
end

% time array
time = (1:length(Grab_raw))/data.streams.(Grab).fs;

% removes the first (t) seconds of signal
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
Grab_raw = Grab_raw(ind:end);
ISOS_raw = ISOS_raw(ind:end);

% Downsample streams and time array by N times%
ISOS_raw = downsample(ISOS_raw, N);
Grab_raw = downsample(Grab_raw, N);
time = downsample(time, N);

% Converts raw mV signal to dFF and detrends to remove photobleaching%
bls = polyfit(ISOS_raw,Grab_raw,1);
Y_fit_all = bls(1) .* ISOS_raw + bls(2);
Y_dF_all = Grab_raw - Y_fit_all; %sF (units mV) is not dFF
Grab_dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(Grab_dFF));
Grab_dFF_detrend = detrend(Grab_dFF);
Grab_dFF_zscore = zscore(Grab_dFF_detrend);

ind1 = find(time > x,1);
ind2 = find(time > y,1);
Grab_snip = Grab_dFF_detrend(ind1:ind2);
time_snip = time(ind1:ind2);
% Creates figure of raw, dFF, and detrended dFF GCaMP signal

f1 = figure;
subplot(4,1,1)
plot(time,Grab_raw,'b')
title(TITLE)
ylabel(sprintf('Raw %s (mV)',Grab_Name),"FontSize",fontSize)
subplot(4,1,2)
plot(time,Grab_dFF,'b')
ylabel(sprintf('%s (dFF)',Grab_Name),"FontSize",fontSize)
subplot(4,1,3)
plot(time,Grab_dFF_detrend,'b')
ylabel(sprintf('%s Detrended (dFF)',Grab_Name),"FontSize",fontSize)
subplot(4,1,4)
plot(time_snip,Grab_snip,'b')
title('10s Window')
ylabel(sprintf('%s dFF',Grab_Name),"FontSize",fontSize)
xlabel('Time (s)')

f2 = figure;
subplot(3,1,1)
yyaxis left
plot(time,Grab_raw,'b')
ylabel(sprintf('Uncorrected %s (mV)',Grab_Name),"FontSize",fontSize)
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca,'YLabel'),'Color','black')
yyaxis right
plot(time,ISOS_raw,'r')
ylabel('Uncorrected Isosbestic (mV)',"FontSize",fontSize)
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca, 'YLabel'), 'Rotation', -90, 'Color','black') % Rotate the right ylabel
legend(sprintf('%s',Grab_Name),'Isosbestic')
title(TITLE,"FontSize",12)
set(gca, 'Children', flipud(get(gca, 'Children')))
subplot(3,1,2)
plot(time,Grab_dFF_detrend,'b')
ylabel(sprintf('%s (dFF)',Grab_Name),"FontSize",fontSize)
subplot(3,1,3)
plot(time,Grab_dFF_zscore,'b')
ylabel(sprintf('%s (z-score)',Grab_Name),"FontSize",fontSize)
xlabel('Time (s)')

% Saves figure as desired file type to user chosen directory
file_name1 = strcat(figsavepath,name,'_fig1',figsavetype);
file_name2 = strcat(figsavepath,name,'_fig2',figsavetype);
print(f1,file_name1,'-dpdf','-bestfit');
print(f2,file_name2,'-dpdf','-bestfit');
