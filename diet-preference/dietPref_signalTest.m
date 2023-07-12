clear all;
close all;
clc;
warning off;
% dietPref_sigTest.m created by Brandon L. Oliver, M.A.
%
% Used to quickly test GRAB sensor flourescence from TDT fiberphotometry tanks.
% Plots the raw, dFF, and z-score dFF on seperate plots and saves a
% figure within the tank folder. Can also downsample if necessary
% by changing the value of 'N' below. To start, 
% This script assumes the following file naming convention:
% IDA_Treatment_IDB_Treatment_Task (example: 12F_V_34M_A_SDWD)
% If only one animal was recorded, then the file should have the following
% naming convention: 12F_V_Empty_NA_SDWD if stream C is empty, or 
% Empty_NA_34M_A_SDWD if stream A is empty.
VERSION = 'v1.0';
%%%%%%%%%%%%%%%%%%%%%%%%% Variables to Change %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRAB_Sensor = 'GrabDA';
figsavetype = '.pdf'; % can change to '.jpg', '.fig', etc.
t = 30; % first t seconds are discarded to remove LED on artifact
N = 10; % downsample signal N times
fontSize = 8; % font size for figure ylabels
experimentLogic = {'VC', 'Lox'}; % Ignores tanks not including these
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("dietPref_signalTest Version: %s\n",VERSION)
myDir = uigetdir(pwd,"Select a folder containing one or more tanks"); 
if myDir == 0
    disp("Select a folder containing one or more tanks")
    return
end
errorFiles = {};
errorMessages = {};
myFiles = dir(myDir);
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
myFiles = myFiles(~endsWith({myFiles.name},{'.jpg','.m','.asv','.xlsx'}));
numFiles = length(myFiles);

for batch = 1:numFiles
    try
    TANK_NAME = fullfile(myDir, myFiles(batch).name);
    figsavepath = strcat(TANK_NAME,'/');
    [~,name,~] = fileparts(TANK_NAME);
    fprintf('Loading %s (%d of %d)...\n',name,batch,numFiles)
    brokenID = strsplit(name,'_');
    emptyLogic = 'Empty';
    emptyIDA = brokenID(1);
    emptyIDC = brokenID(3);
    if ~strcmp(brokenID(2),experimentLogic(1)) && ~strcmp(brokenID(2),experimentLogic(2))
        disp('File is not a part of current experiment. Skipping.')
        continue
    end
    for i = 1:2
        if i == 1
            if strcmp(emptyIDA,emptyLogic)
                disp('No data in stream A.')
                continue
            end
            ISOS = 'x405A';
            GRAB = 'x465A';
            ID = brokenID(1);
            treatment = brokenID(2);
            task = brokenID(5);
            TITLE = strcat(ID,{' '},treatment,{' '},task);
            data = TDTbin2mat(TANK_NAME, 'TYPE', {'streams'});
            ISOS_raw = data.streams.(ISOS).data;
            GRAB_raw = data.streams.(GRAB).data;
        end
        if i == 2
           if strcmp(emptyIDC,emptyLogic)
                disp('No data in stream C.')
                continue
            end
            ISOS = 'x405C';
            GRAB = 'x465C';
            ID = brokenID(3);
            treatment = brokenID(4);
            task = brokenID(5);
            TITLE = strcat(ID,{' '},treatment,{' '},task);
            data = TDTbin2mat(TANK_NAME, 'TYPE', {'streams'});
            ISOS_raw = data.streams.(ISOS).data;
            GRAB_raw = data.streams.(GRAB).data; 
        end
        % Checks for unequal isosbestic and GRAB sensor signal length correcting if
        % necessary
        if length(GRAB_raw) < length(ISOS_raw)
            disp('Isosbestic signal array is longer than GRAB signal array')
            ISOS_raw = ISOS_raw(1:length(GRAB_raw));
            disp('Corrected.')
        elseif length(GRAB_raw) > length(ISOS_raw)
            disp('Isosbestic signal array is shorter than GRAB signal array')
            GRAB_raw = GRAB_raw(1:length(ISOS_raw));
            disp('Corrected.')
        else
            disp('Isosbestic and GRAB signal arrays are of equal size')
            disp('No correction necessary.')
        end
        
        % time array
        time = (1:length(GRAB_raw))/data.streams.(GRAB).fs;
        
        % removes the first (t) seconds of signal
        ind = find(time>t,1);% find first index of when time crosses threshold
        time = time(ind:end); % reformat vector to only include allowed time
        GRAB_raw = GRAB_raw(ind:end);
        ISOS_raw = ISOS_raw(ind:end);
        
        % Downsample streams and time array by N times%
        ISOS_raw = downsample(ISOS_raw, N);
        GRAB_raw = downsample(GRAB_raw, N);
        time = downsample(time, N);
        
        % Converts raw mV signal to dFF and detrends to remove photobleaching%
        bls = polyfit(ISOS_raw,GRAB_raw,1);
        Y_fit_all = bls(1) .* ISOS_raw + bls(2);
        Y_dF_all = GRAB_raw - Y_fit_all; %sF (units mV) is not dFF
        GRAB_dFF = 100*(Y_dF_all)./Y_fit_all;
        std_dFF = std(double(GRAB_dFF));
        GRAB_dFF = detrend(GRAB_dFF);
        GRAB_dFF_z = zscore(GRAB_dFF);

        % Calculates dFF and zscore for isosbestic
        ISOS_raw = detrend(ISOS_raw);
        ISOS_mean = mean(ISOS_raw);
        ISOS_dF = ISOS_raw - ISOS_mean;
        ISOS_dFF = ISOS_dF ./ ISOS_mean;
        ISOS_dFF_z = zscore(ISOS_dFF);
        
        f1 = figure('Visible','off');
        
        subplot(3,1,1)
        title(TITLE,"FontSize",12)
        yyaxis left
        plot(time,GRAB_raw,'b')
        ylabel(sprintf('%s (mV)',GRAB_Sensor),"FontSize",fontSize)
        set(gca,'YColor','black','Box', 'on','Color','w')
        set(get(gca,'YLabel'),'Color','black')
        xlim([t floor(time(1,end))])
        ylim1 = min(GRAB_raw)*(0.99);
        ylim2 = max(GRAB_raw)*(1.05);
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
        f1a = legend(sprintf('%s (465nm)',GRAB_Sensor), 'Isosbestic (405nm)','Orientation','horizontal','Location','best');
        f1a.FontSize = fontSize;
        
        
        subplot(3,1,2)
        yyaxis left
        plot(time,GRAB_dFF,'b')
        ylabel([sprintf('%s\n',GRAB_Sensor),' \DeltaF/F_{0}'],"FontSize",fontSize,'Interpreter','tex')
        set(gca,'YColor','black','Box', 'on','Color','w')
        set(get(gca,'YLabel'),'Color','black')
        xlim([t floor(time(1,end))])
        ylim1 = min(GRAB_dFF)*(2.1);
        ylim2 = max(GRAB_dFF)*(1.1);
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
        f1b = legend(sprintf('%s (465nm)',GRAB_Sensor), 'Isosbestic (405nm)','Orientation','horizontal','Location','best');
        f1b.FontSize = fontSize;
        
        subplot(3,1,3)
        yyaxis left
        plot(time,GRAB_dFF_z,'b')
        xlim([t floor(time(1,end))])
        ylim1 = min(GRAB_dFF_z)*(2.1);
        ylim2 = max(GRAB_dFF_z)*(1.1);
        ylim([ylim1,ylim2])
        ylabel([sprintf('%s\n',GRAB_Sensor),' \DeltaF/F_{0} (z-Score)'],'FontSize',fontSize,'Interpreter','tex')
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
        f1c = legend(sprintf('%s (465nm)',GRAB_Sensor), 'Isosbestic (405nm)','Orientation','horizontal','Location','best');
        f1c.FontSize = fontSize;
        
        
        
        file_name1 = char(strcat(figsavepath,TITLE,figsavetype));
        print(f1,file_name1,'-dpdf','-bestfit');
        close(f1)
        disp('Done.')
    end
    catch ex
        fprintf('Encountered an error with file %d\n',name)
        % Store the file name and error message
        errorFiles{end+1} = name;
        errorMessages{end+1} = ex.message;
    end
end
close all
for j = 1:numel(errorFiles)
    if numel(errorFiles) == 0
        disp('No file errors occured!')
        break
    end
    fprintf('Error in file: %s\n', errorFiles{i});
    fprintf('Error message: %s\n', errorMessages{i});
end