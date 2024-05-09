clear all;
close all;
clc;
warning off;
VERSION = 'v2.0';
% Grab_signalTest.m created by Brandon L. Oliver, M.A.
% Used to quickly test Grab sensor flourescence from TDT fiberphotometry tanks.
% Plots the raw, dFF, and detrended dFF on seperate plots and saves a
% figure to a user designated directory. Can also downsample if necessary
% by changing the value of 'N' below.
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figsavetype = '.pdf'; % can change to '.jpg', '.fig', etc.
t = 10; % first t seconds are discarded to remove LED on artifact
N = 1000; % downsample signal N times
channel = 2; % 1 = A, 2 = C
session_durration = 3600;
fontSize = 8; % font size for figure ylabels
figureSize = [100,100,1500,800]; % Set the desired figure size
figSnip = [2000, 2060];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Grab_signalTest %s\n',VERSION)
% Gets tank from UI pop-up window
TANK_NAME = uigetdir('/Volumes/KINGSTON/FIBRE PHOTOMETRY - SYNAPSE', 'Select a tank to plot');
if TANK_NAME == 0
    disp('Select a file to start!')
    return
end

figsavepath = strcat(TANK_NAME,'/');
[~,name,~] = fileparts(TANK_NAME);
brokenID = strsplit(name,'_');
if channel == 1
    ISOS = 'x405A'; % set name of isosbestic signal
    Grab = 'x465A'; % set name of Grab signal
    ROI = 'DLS';
    ID = brokenID(2);
    Grab_Sensor = 'GrabDA'; % example: 'GrabDA4.4'
elseif channel == 2
    ISOS = 'x405C'; % set name of isosbestic signal
    Grab = 'x465C'; % set name of Grab signal
    ROI = 'NAc';
    ID = brokenID(2);
    Grab_Sensor = 'GrabDA'; % example: 'GrabDA4.4'
else
    disp('Cannot find isosbestic signal. Check the naming and try again.')
end


TITLE = strcat(ID,{' '},Grab_Sensor,{' '},ROI);
data = TDTbin2mat(TANK_NAME, 'T2', session_durration,'TYPE', {'streams','epocs'});
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
Grab_dFF = detrend(Grab_dFF);
Grab_dFF_z = zscore(Grab_dFF);


% Converts raw mV isosbestic signal to dFF and zScore for plotting %
ISOS_dFF = detrend(deltaFF(ISOS_raw));
ISOS_dFF_z = zscore(ISOS_dFF);



f1 = figure;
set(f1,'Position',figureSize);

subplot(3,1,1)
title(TITLE,"FontSize",12)
yyaxis left
plot(time,Grab_raw,'b')
ylabel(sprintf('%s (mV)',Grab_Sensor),"FontSize",fontSize)
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca,'YLabel'),'Color','black')
xlim([t floor(time(1,end))])
ylim1 = min(Grab_raw)*(0.99);
ylim2 = max(Grab_raw)*(1.05);
ylim([ylim1,ylim2])

yyaxis right
plot(time,ISOS_raw,'r')
ylabel(sprintf('Isosbestic (mV)'),"FontSize",fontSize)
xlim([t floor(time(1,end))])
ylim1 = min(ISOS_raw)*(0.95);
ylim2 = max(ISOS_raw)*(1.01);
ylim([ylim1,ylim2])
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca, 'YLabel'), 'Rotation', -90, 'Color','black') % Rotate the right ylabel
f1a = legend(sprintf('%s (465nm)',Grab_Sensor), 'Isosbestic (405nm)','Orientation','horizontal','Location','northeast');
f1a.FontSize = fontSize;


subplot(3,1,2)
yyaxis left
plot(time,Grab_dFF,'b')
ylabel([sprintf('%s\n',Grab_Sensor),' \DeltaF/F_{0}'],"FontSize",fontSize,'Interpreter','tex')
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca,'YLabel'),'Color','black')
xlim([t floor(time(1,end))])
ylim1 = min(Grab_dFF)*(2.1);
ylim2 = max(Grab_dFF)*(1.1);
ylim([ylim1,ylim2])

yyaxis right
plot(time,ISOS_dFF,'r')
ylabel('Isosbestic \DeltaF/F_{0}','FontSize',fontSize,'Interpreter','tex')
xlim([t floor(time(1,end))])
ylim1 = min(ISOS_dFF)*(1.1);
ylim2 = max(ISOS_dFF)*(2.1);
ylim([ylim1,ylim2])
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca, 'YLabel'), 'Rotation', -90, 'Color','black') % Rotate the right ylabel
f1b = legend(sprintf('%s (465nm)',Grab_Sensor), 'Isosbestic (405nm)','Orientation','horizontal','Location','northeast');
f1b.FontSize = fontSize;

subplot(3,1,3)
yyaxis left
plot(time,Grab_dFF_z,'b')
xlim([t floor(time(1,end))])
ylim1 = min(Grab_dFF_z)*(2.1);
ylim2 = max(Grab_dFF_z)*(1.1);
ylim([ylim1,ylim2])
ylabel([sprintf('%s\n',Grab_Sensor),' \DeltaF/F_{0} (z-Score)'],'FontSize',fontSize,'Interpreter','tex')
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca,'YLabel'),'Color','black')
xlabel('Time (s)')

yyaxis right
plot(time,ISOS_dFF_z,'r')
ylabel('Isosbestic \DeltaF/F_{0} (z-Score)','FontSize',fontSize,'Interpreter','tex')
xlim([t floor(time(1,end))])
ylim1 = min(ISOS_dFF_z)*(1.1);
ylim2 = max(ISOS_dFF_z)*(2.1);
ylim([ylim1,ylim2])
set(gca,'YColor','black','Box', 'on','Color','w')
set(get(gca, 'YLabel'), 'Rotation', -90, 'Color','black') % Rotate the right ylabel
f1c = legend(sprintf('%s (465nm)',Grab_Sensor), 'Isosbestic (405nm)','Orientation','horizontal','Location','northeast');
f1c.FontSize = fontSize;



file_name1 = char(strcat(figsavepath,TITLE,' fig1',figsavetype));
print(f1,file_name1,'-dpdf','-vector','-bestfit');

% peakDist = 0.2;
% ind = find(time>data.epocs.coc_.onset(1,1),1);
% ind2 = find(time>120,1);
% ind2end = find(time>240,1);
% ind3 = find(time>420,1);
% ind3end = find(time>540,1);
% 
% 
% MAD1 = mad(Grab_dFF_z(1,1:ind), 1);
% MAD2 = mad(Grab_dFF_z(1,ind:end),1);
% [pks,locs,w,p] = findpeaks(Grab_dFF_z(1,1:ind), time(1,1:ind), 'MinPeakHeight', MAD1,...
%             "MinPeakDistance", peakDist);
% [pks2,locs2,w2,p2] = findpeaks(Grab_dFF_z(1,ind:end), time(1,ind:end), 'MinPeakHeight', MAD2,...
%             "MinPeakDistance", peakDist);
% 
% b4cocpks = length(pks);
% b4cocpks_m = (b4cocpks/time(1,end))*60;
% aftrcocpks = length(pks2);
% aftrcocpks_m = (aftrcocpks/time(1,end))*60;
% 
% pkanalysis = table(b4cocpks,aftrcocpks,b4cocpks_m,aftrcocpks_m,'VariableNames',{'NumPeaks Before Coc', 'NumPeaks After Coc', ...
%     'Peaks/m Before Coc', 'Peaks/m After Coc'});
% 
% f2 = figure;
% subplot(4,1,1)
% plot(time,ISOS_raw, 'r')
% xlim([t floor(time(1,end))])
% ylabel('Isosbestic Uncorrected (mV)','FontSize',fontSize)
% subplot(4,1,2)
% plot(time,Grab_raw,'b')
% xlim([t floor(time(1,end))])
% ylabel('GrabDA Uncorrected (mV)','FontSize',fontSize)
% subplot(4,1,3)
% plot(time,ISOS_dFF_z,'r')
% xlim([t floor(time(1,end))])
% ylabel('Isosbestic \DeltaF/F_{0} (z-Score)','FontSize',fontSize,'Interpreter','tex')
% subplot(4,1,4)
% plot(time,Grab_dFF_z,'b')
% title(sprintf('Peaks/m Before Cocaine: %.2f, After Cocaine: %.2f',b4cocpks_m,aftrcocpks_m))
% xline(data.epocs.coc_.onset(1,1),'-red', 'Infusion','LineWidth',2,'FontSize',16)
% xlim([t floor(time(1,end))])
% ylabel('GrabDA \DeltaF/F_{0} (z-Score)','FontSize',fontSize,'Interpreter','tex')
% xlabel('Time (s)')
% 
% f3 = figure;
% subplot(2,1,1)
% plot(time(1,ind2:ind2end),Grab_dFF_z(1,ind2:ind2end))
% title('Before Cocaine (2 mins)')
% xlim([floor(time(1,ind2)) floor(time(1,ind2end))])
% ylim([-4,6])
% ylabel('GrabDA \DeltaF/F_{0} (z-Score)','FontSize',fontSize,'Interpreter','tex')
% xlabel('Time (s)')
% subplot(2,1,2)
% plot(time(1,ind3:ind3end), Grab_dFF_z(1,ind3:ind3end))
% title('After Cocaine (2 min)')
% xlim([floor(time(1,ind3)) floor(time(1,ind3end))])
% ylim([-4,6])
% ylabel('GrabDA \DeltaF/F_{0} (z-Score)','FontSize',fontSize,'Interpreter','tex')
% xlabel('Time (s)')

idx = find(time>figSnip(1),1);
idx2 = find(time>figSnip(2),1);
timeSnip = time(1, idx:idx2);
grabSnip = Grab_dFF_z(1, idx:idx2);
ts1 = 0 + (1:length(grabSnip)) / data.streams.(Grab).fs*N;
f4 = figure;
set(f4,'Position',[800,550,1200,420]);
plot(timeSnip,grabSnip, 'Color', 'b')
ylabel('GrabDA \DeltaF/F_{0} (z-Score)','Interpreter','tex')
xlabel('Time (s)')
title(TITLE)
subtitle(sprintf('%d Second Window', (figSnip(2)-figSnip(1))))
xlim([figSnip(1),figSnip(2)])
disp('Done.')
