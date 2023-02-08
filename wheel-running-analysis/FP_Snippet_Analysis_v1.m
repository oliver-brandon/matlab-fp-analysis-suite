%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                            
%FP_Snippet_Analysis.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)

VERSION = "1.0";

%House Keeping%
clear all; clc; close all;
%removes polyfit warning%
warning('off','all');

%Choose batch or single tank analysis%
tanks2analyze = 1;%1 = batch, 2 = single%
%Choose zscore paramater%
%1 = convert stream to zscore before snipping, 2 = convert snippets to
%zscore (local zscore), 3 = choose your own zscore window
snip_z = 1;
z_window = [660,5400];
%Choose minimum peak distance%
peakDist = 0.2;
%Choose session duration%
session_duration = 3600; %duration of recording in seconds
%amount to downsample signals%
N = 100; %Downsample N times
%Time window for snippets%
%Currently, the code needs exactly three snippet windows to run%
snip1 = [500,800];
snip2 = [750,1050];
snip3 = [1000,1300];
% snip4 = [3000,3300];
% snip5 = [3000,3300];
snip_list = {snip1,snip2,snip3};
%seconds before snip to calculate MAD from%
usemadwindow = 3;%1=use window to calculate MAD, 2=calculate MAD for whole session,
%3=use MAD values from another session (save MAD values as a cell .mat)
DLS_PostVehMAD = "DLS_Post_Veh_U_MAD.mat";
NAc_PostVehMAD = "NAC_Post_Veh_U_MAD.mat";
madtime = [660,1740];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(VERSION)
if tanks2analyze == 1
    myDir = uigetdir; %gets directory%
    myFiles = dir(myDir); %gets all tanks in directory%
    myFiles = myFiles(~ismember({myFiles.name},{'.','..'}));
    numFiles = length(myFiles);
    disp("Starting batch analysis...")
    for a = 1:numel(snip_list)
        snip = snip_list{a};
        fprintf('Snipping streams from %d to %d seconds...\n',snip(:,1),snip(:,2))

        for i = 1:length(myFiles)
            BLOCKPATH = fullfile(myDir, myFiles(i).name);
            [tankpath,tankname,ext] = fileparts(BLOCKPATH);
            fprintf('Loading tank %d of %d...\n',i,length(myFiles))
            data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});
            
            %Stream Stores%
            DLS_ISOS = 'x405A'; % name of the 405A store
            DLS_DA = 'x465A'; % name of the 465A store
            NAc_ISOS = 'x405C'; % name of the 405C store
            NAc_DA = 'x465C'; % name of the 465C store
            
         
            %time array used for all streams%
            time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
            %removes the first (t) seconds where the data is wild due to turning on LEDs%
            t = 20; % time threshold below which we will discard
            ind = find(time>t,1);% find first index of when time crosses threshold
            time = time(ind:end); % reformat vector to only include allowed time
            data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
            data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
            data.streams.(NAc_DA).data = data.streams.(NAc_DA).data(ind:end);
            data.streams.(NAc_ISOS).data = data.streams.(NAc_ISOS).data(ind:end);
            
            %downsample streams and time array by N times%
            data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
            data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
            data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
            data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
            time = downsample(time, N);
            
            %detrend & dFF%
            %465A%
            bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
            Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
            Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
            dFF = 100*(Y_dF_all)./Y_fit_all;
            std_dFF = std(double(dFF));
           
            
            %465C%
            bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
            Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
            Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
            dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
            std_dFF2 = std(double(dFF2));
            
            
            %%%converts streams to zscore before snipping them%%%
            if snip_z == 1
                detrend_465A = detrend(dFF);
                detrend_465A = zscore(detrend_465A);
%                 detrend_465A = smoothdata(detrend_465A,'lowess');
                detrend_465C = detrend(dFF2);
                detrend_465C = zscore(detrend_465C);
%                 detrend_465C = smoothdata(detrend_465C,'lowess');
            elseif snip_z == 2
                detrend_465A = detrend(dFF);
%                 detrend_465A = smoothdata(detrend_465A,'lowess');
                detrend_465C = detrend(dFF2);
%                 detrend_465C = smoothdata(detrend_465C,'lowess');
            elseif snip_z == 3
                detrend_465A = detrend(dFF);
                detrend_465C = detrend(dFF2);
                z_win_start = find(time>z_window(:,1),1);
                z_win_end = find(time>z_window(:,2),1);
                detrend_465A = zscore(detrend_465A(z_win_start:z_win_end));
                detrend_465C = zscore(detrend_465C(z_win_start:z_win_end));
            end
           
            %calculates and plots median absolute deviation for both 465 signals%
            MAD1 = mad(detrend_465A, 1);
            MAD2 = mad(detrend_465C, 1);
            DLS_MAD_STORE{i} = MAD1;
            NAC_MAD_STORE{i} = MAD2;
            
            [pks,locs,w,p] = findpeaks(detrend_465A, time, 'MinPeakProminence', MAD1,...
                "MinPeakDistance", peakDist);
            [pks2,locs2,w2,p2] = findpeaks(detrend_465C, time, 'MinPeakProminence', MAD2,...
                "MinPeakDistance", peakDist);
            
            DLS_pks = length(pks);
            DLS_pk_min = (DLS_pks/session_duration)*60;
            DLS_amp_max = max(pks);
            DLS_amp_avg = mean(pks);
            NAc_pks = length(pks2);
            NAc_pk_min = (NAc_pks/session_duration)*60;
            NAc_amp_max = max(pks2);
            NAc_amp_avg = mean(pks2);
            
            
            snipstart = find(time>snip(:,1),1);
            snipend = find(time>snip(:,2),1);
            
           
            DLS_trim = detrend_465A(snipstart:snipend);
            NAc_trim = detrend_465C(snipstart:snipend);
            time_trim = time(snipstart:snipend);
            if snip_z == 1
                disp('')
            elseif snip_z == 2
               %converts snippets to zscore%
                DLS_trim = zscore(DLS_trim);
                NAc_trim = zscore(NAc_trim);
            elseif snip_z == 3
                disp('')
            end
            
            
            if usemadwindow == 1
                %for MAD window%
                madstart = find(time>madtime(:,1),1);
                madend = find(time>madtime(:,2),1);
                MADwindow1 = mad(zscore(detrend_465A(madstart:madend)),1);
                MADwindow2 = mad(zscore(detrend_465C(madstart:madend)),1);
                [snip_pks,snip_locs] = findpeaks(DLS_trim,time_trim,'MinPeakProminence',MADwindow1,...
                        "MinPeakDistance", peakDist);
                [snip_pks2,snip_locs2] = findpeaks(NAc_trim,time_trim,'MinPeakProminence',MADwindow2,...
                        "MinPeakDistance", peakDist);
            elseif usemadwindow == 2
                [snip_pks,snip_locs] = findpeaks(DLS_trim,time_trim,'MinPeakProminence',MAD1,...
                    "MinPeakDistance", peakDist);
                [snip_pks2,snip_locs2] = findpeaks(NAc_trim,time_trim,'MinPeakProminence',MAD2,...
                    "MinPeakDistance", peakDist);
            elseif usemadwindow ==3
                DLS_MAD_VEH = load("DLS_Post_Veh_U_MAD.mat");
                NAC_MAD_VEH = load("NAC_Post_Veh_U_MAD.mat");
                DLS_MAD_VEH = DLS_MAD_VEH.DLS_MAD_STORE;
                NAC_MAD_VEH = NAC_MAD_VEH.NAC_MAD_STORE;
                [snip_pks,snip_locs] = findpeaks(DLS_trim,time_trim,'MinPeakProminence',DLS_MAD_VEH{i},...
                    "MinPeakDistance", peakDist);
                [snip_pks2,snip_locs2] = findpeaks(NAc_trim,time_trim,'MinPeakProminence',NAC_MAD_VEH{i},...
                    "MinPeakDistance", peakDist);
            end
            

            DLS_trim_pks(i,1) = length(snip_pks);
            DLS_trim_freq(i,1) = (DLS_trim_pks(i,1)/(snip(:,2)-snip(:,1)))*60;
            NAc_trim_pks(i,1) = length(snip_pks2);
            NAc_trim_freq(i,1) = (NAc_trim_pks(i,1)/(snip(:,2)-snip(:,1)))*60;
            DLS_trim_amp(i,1) = mean(snip_pks);
            NAc_trim_amp(i,1) = mean(snip_pks2);
            DLS_stream_store(i,:) = detrend_465A;
            NAc_stream_store(i,:) = detrend_465C;
            
            if a == 1
                snip_peak_analysis = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
                DLS_snip_store(i,:) = DLS_trim;
                NAc_snip_store(i,:) = NAc_trim;
            elseif a == 2
                snip_peak_analysis2 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
                DLS_snip_store2(i,:) = DLS_trim;
                NAc_snip_store2(i,:) = NAc_trim;
            elseif a == 3
                snip_peak_analysis3 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
                DLS_snip_store3(i,:) = DLS_trim;
                NAc_snip_store3(i,:) = NAc_trim;
            elseif a == 4
            snip_peak_analysis4 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
            DLS_snip_store4(i,:) = DLS_trim;
            NAc_snip_store4(i,:) = NAc_trim;
            elseif a == 5
            snip_peak_analysis5 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
            DLS_snip_store5(i,:) = DLS_trim;
            NAc_snip_store5(i,:) = NAc_trim;
            end

        end

    end
    
    elseif tanks2analyze == 2
        myFiles = uigetdir;
        BLOCKPATH = myFiles;
        numFiles = 1;
        disp("Starting single tank analysis...")
        for a = 1:numel(snip_list)
            snip = snip_list{a};
            fprintf('Snipping streams from %d to %d seconds...\n',snip(:,1),snip(:,2))
            data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});
            
            %Stream Stores%
            DLS_ISOS = 'x405A'; % name of the 405A store
            DLS_DA = 'x465A'; % name of the 465A store
            NAc_ISOS = 'x405C'; % name of the 405C store
            NAc_DA = 'x465C'; % name of the 465C store
            N = 100; %Downsample N times
            
     
            %time array used for all streams%
            time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
            %removes the first (t) seconds where the data is wild due to turning on LEDs%
            t = 20; % time threshold below which we will discard
            ind = find(time>t,1);% find first index of when time crosses threshold
            time = time(ind:end); % reformat vector to only include allowed time
            data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
            data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
            data.streams.(NAc_DA).data = data.streams.(NAc_DA).data(ind:end);
            data.streams.(NAc_ISOS).data = data.streams.(NAc_ISOS).data(ind:end);
            
            %downsample streams and time array by N times%
            data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
            data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
            data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
            data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
            time = downsample(time, N);
            
            %detrend & dFF%
            %465A%
            bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
            Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
            Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
            dFF = 100*(Y_dF_all)./Y_fit_all;
            std_dFF = std(double(dFF));
           
            
            %465C%
            bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
            Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
            Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
            dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
            std_dFF2 = std(double(dFF2));
            
            
            %%%converts streams to zscore before snipping them%%%
            if snip_z == 1
                detrend_465A = detrend(dFF);
                detrend_465A = zscore(detrend_465A);
    %                 detrend_465A = smoothdata(detrend_465A,'lowess');
                detrend_465C = detrend(dFF2);
                detrend_465C = zscore(detrend_465C);
    %                 detrend_465C = smoothdata(detrend_465C,'lowess');
            elseif snip_z == 2
                detrend_465A = detrend(dFF);
    %                 detrend_465A = smoothdata(detrend_465A,'lowess');
                detrend_465C = detrend(dFF2);
    %                 detrend_465C = smoothdata(detrend_465C,'lowess');
            elseif snip_z == 3
                detrend_465A = detrend(dFF);
                detrend_465C = detrend(dFF2);
                z_win_start = find(time>z_window(:,1),1);
                z_win_end = find(time>z_window(:,2),1);
                detrend_465A = zscore(detrend_465A(z_win_start:z_win_end));
                detrend_465C = zscore(detrend_465C(z_win_start:z_win_end));
            end
           
            %calculates and plots median absolute deviation for both 465 signals%
            MAD1 = mad(detrend_465A, 1);
            MAD2 = mad(detrend_465C, 1);
            
            [pks,locs,w,p] = findpeaks(detrend_465A, time, 'MinPeakProminence', MAD1,...
                "MinPeakDistance", peakDist);
            [pks2,locs2,w2,p2] = findpeaks(detrend_465C, time, 'MinPeakProminence', MAD2,...
                "MinPeakDistance", peakDist);
            
            DLS_pks = length(pks);
            DLS_pk_min = (DLS_pks/session_duration)*60;
            DLS_amp_max = max(pks);
            DLS_amp_avg = mean(pks);
            NAc_pks = length(pks2);
            NAc_pk_min = (NAc_pks/session_duration)*60;
            NAc_amp_max = max(pks2);
            NAc_amp_avg = mean(pks2);
            
            
            snipstart = find(time>snip(:,1),1);
            snipend = find(time>snip(:,2),1);
            
           
            DLS_trim = detrend_465A(snipstart:snipend);
            NAc_trim = detrend_465C(snipstart:snipend);
            time_trim = time(snipstart:snipend);
            if snip_z == 1
                disp('')
            elseif snip_z == 2
               %converts snippets to zscore%
                DLS_trim = zscore(DLS_trim);
                NAc_trim = zscore(NAc_trim);
            elseif snip_z == 3
                disp('')
            end
            
            
            if usemadwindow == 1
                %for MAD window%
                madstart = find(time>madtime(:,1),1);
                madend = find(time>madtime(:,2),1);
                MADwindow1 = mad(zscore(detrend_465A(madstart:madend)),1);
                MADwindow2 = mad(zscore(detrend_465C(madstart:madend)),1);
                [snip_pks,snip_locs] = findpeaks(DLS_trim,time_trim,'MinPeakProminence',MADwindow1,...
                        "MinPeakDistance", peakDist);
                [snip_pks2,snip_locs2] = findpeaks(NAc_trim,time_trim,'MinPeakProminence',MADwindow2,...
                        "MinPeakDistance", peakDist);
            elseif usemadwindow == 2
                [snip_pks,snip_locs] = findpeaks(DLS_trim,time_trim,'MinPeakProminence',MAD1,...
                    "MinPeakDistance", peakDist);
                [snip_pks2,snip_locs2] = findpeaks(NAc_trim,time_trim,'MinPeakProminence',MAD2,...
                    "MinPeakDistance", peakDist);
            end
            
    
            DLS_trim_pks = length(snip_pks);
            DLS_trim_freq = (DLS_trim_pks/(snip(:,2)-snip(:,1)))*60;
            NAc_trim_pks = length(snip_pks2);
            NAc_trim_freq = (NAc_trim_pks/(snip(:,2)-snip(:,1)))*60;
            DLS_trim_amp = mean(snip_pks);
            NAc_trim_amp = mean(snip_pks2);
            DLS_stream_store = detrend_465A;
            NAc_stream_store = detrend_465C;
            
            if a == 1
                snip_peak_analysis = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
                DLS_snip_store = DLS_trim;
                NAc_snip_store = NAc_trim;
            elseif a == 2
                snip_peak_analysis2 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
                DLS_snip_store2 = DLS_trim;
                NAc_snip_store2 = NAc_trim;
            elseif a == 3
                snip_peak_analysis3 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
                DLS_snip_store3 = DLS_trim;
                NAc_snip_store3 = NAc_trim;
            elseif a == 4
            snip_peak_analysis4 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
            DLS_snip_store4 = DLS_trim;
            NAc_snip_store4 = NAc_trim;
            elseif a == 5
            snip_peak_analysis5 = table(DLS_trim_pks,NAc_trim_pks,DLS_trim_amp,NAc_trim_amp, ...
                        DLS_trim_freq,NAc_trim_freq,'VariableNames', {'DLS Peaks','NAc Peaks',...
                        'DLS Avg Amp','NAc Avg Amp','DLS Peaks/m','NAc Peaks/m'});
            DLS_snip_store5 = DLS_trim;
            NAc_snip_store5 = NAc_trim;
            end
            

        end
        disp("Plotting signals...")
        figure;
        subplot(2,1,1)
        title("DLS 465A")
        xlabel("Z-Score")
        plot(time,detrend_465A,'r');
        subplot(2,1,2)
        title("NAc 465C")
        xlabel("Z-Score")
        ylabel("Time (s)")
        plot(time,detrend_465C,'b');
end
    %UITable (figure) that displays "snip_peak_analysis" table%
    figure;
    uitable('Data',snip_peak_analysis{:,:},'ColumnName',snip_peak_analysis.Properties.VariableNames,...
        'RowName',snip_peak_analysis.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    figure;
    uitable('Data',snip_peak_analysis2{:,:},'ColumnName',snip_peak_analysis2.Properties.VariableNames,...
        'RowName',snip_peak_analysis2.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    figure;
    uitable('Data',snip_peak_analysis3{:,:},'ColumnName',snip_peak_analysis3.Properties.VariableNames,...
        'RowName',snip_peak_analysis3.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%     figure;
%     uitable('Data',snip_peak_analysis4{:,:},'ColumnName',snip_peak_analysis4.Properties.VariableNames,...
%         'RowName',snip_peak_analysis4.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
%     figure;
%     uitable('Data',snip_peak_analysis5{:,:},'ColumnName',snip_peak_analysis5.Properties.VariableNames,...
%         'RowName',snip_peak_analysis5.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


% 
% figure;
% subplot(2,1,1)
% plot(time,detrend_465A);
% subplot(2,1,2)
% plot(time,detrend_465C);






% figure
% bar([DLS_trim_pks,NAc_trim_pks])
% loco_peak_analysis = table(DLS_pks, DLS_pk_min, DLS_amp_max, DLS_amp_avg, ...
% NAc_pks, NAc_pk_min, NAc_amp_max, NAc_amp_avg, ...
% 'VariableNames', {'DLS Peaks','DLS Peaks/m','DLS Max Amp', ...
% 'DLS Avg Amp','NAc Peaks','NAc Peaks/m','Nac Max Amp', ...
% 'NAc Avg Amp'});
% 
% %UITable (figure) that displays "peak_analysis" table%
% figure;
% uitable('Data',loco_peak_analysis{:,:},'ColumnName',loco_peak_analysis.Properties.VariableNames,...
%     'RowName',loco_peak_analysis.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
disp("Analysis complete!")
fprintf('Successfully analyzed %d tank(s)!',numFiles)


