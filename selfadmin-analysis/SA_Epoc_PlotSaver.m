%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%SA_Epoc_PlotSaver.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)
clear; clc; close all;
VERSION = "1.1";

%%%% how to use %%%%
% set the variable below named "figure_savepath" to your desired figure
% savepath and the "epoc" variable to the desired TTL for both streams 
% (even if one of the streams is empty). Run the script and a UI window 
% will appear asking you to choose a folder containing one or more TDT 
% fiber photometry tanks (if plotting one tank, you still need to place 
% the tank in an empty folder). The tanks need to have the following naming
% convention: IDA_Task-Day_IDB_TaskDay

%%%% returns %%%%
% Saves images of the average signal (+/- SEM) aligned to an epoc of choice
% Note: the default is .jpg, but you can change this by changing the 
% variable "savetype" below to your desired filetype.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% edit these variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure_savepath = '/Users/brandon/personal-drive/self_admin/coc_sa/PrL-aIC/FR1-Cocaine/figs/'; % include forward slash at end of path
% epoc = {'aRw/','bRw/'};
% epoc = {'aRL/','bRL/'};
% epoc = {'aHL/','bHL/'};
epoc = {'aReward','bReward'};
% epoc = {'aActiveRew','bActiveRew'};
% epoc = {'aActiveTimeout','bActiveTimeout'};
% epoc = {'aRewTimeout','bRewTimeout'};
EPOCNAME = 'Infusion';
savetype = ".pdf";
TRANGE = [-2 7]; %window size [start time relative to epoc onset, entire duration]
BASELINE_PER = [-3 -1]; % baseline period before stim
N = 10; % Downsample Nx
addVertLine = 0; % 1 = yes, 0 = no
vertLineX = 2; % if adding vertical line, specifies the x coord
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDir = uigetdir('/Users/brandon/personal-drive/self_admin/coc_sa/PrL-aIC/FR1-Cocaine/tanks',"Select a folder containing one or more tanks"); 
fprintf("SA_Epoc_PlotSaver Version: %s\n",VERSION)
tic
if myDir == 0
    disp("Select a folder containing one or more tanks")
    return
end
myFiles = dir(myDir);
myFiles = myFiles(~startsWith({myFiles.name},{'.','..','_'}));
myFiles = myFiles(~endsWith({myFiles.name},{'.jpg','.m','.asv'}));
fig = cell(length(myFiles),1);
TITLE = cell(length(myFiles),1);
if length(myFiles) > 1
    disp("Starting batch tank plot...")
elseif length(myFiles) == 1
    disp("Starting single tank plot")
end

for batch = 1:length(myFiles)
    fprintf('Loading tank %d of %d...\n',batch,length(myFiles))
    BLOCKPATH = fullfile(myDir, myFiles(batch).name);
    [~,name,~] = fileparts(BLOCKPATH);
    emptyID = 'Empty';
    brokenID = strsplit(name,'_');
    animalIDA = char(brokenID{1});
    animalIDC = char(brokenID{3});
    taskA = char(brokenID{2});
    taskC = char(brokenID{4});
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs','streams'});
    for streamAorC = 1:2
        if streamAorC == 1
            emptylogicA = strcmp(animalIDA,emptyID);
            if emptylogicA == 1
                disp("Stream A is empty")
                continue
            elseif emptylogicA == 0
                [rewardTimestamps, rewardTimeout, timeoutTimestamps] = separateActivePoke(data.epocs.aRL_.onset, 10);
                [data] = createEpoc(data, rewardTimestamps, 'aActiveRew');
                [data] = createEpoc(data, rewardTimeout, 'aRewTimeout');
                [data] = createEpoc(data, timeoutTimestamps, 'aActiveTimeout');
                [data] = createEpoc(data, data.epocs.aRw_.offset, 'aReward');

            end
        end
        if streamAorC == 2
            emptylogicC = strcmp(animalIDC,emptyID);
            if emptylogicC == 1
                disp("Stream C is empty")
                continue
            elseif emptylogicC == 0
                [rewardTimestamps, rewardTimeout, timeoutTimestamps] = separateActivePoke(data.epocs.bRL_.onset, 10);
                [data] = createEpoc(data, rewardTimestamps, 'bActiveRew');
                [data] = createEpoc(data, rewardTimeout, 'bRewTimeout');
                [data] = createEpoc(data, timeoutTimestamps, 'bActiveTimeout');
                [data] = createEpoc(data, data.epocs.bRw_.offset, 'bReward');
            end
        end
        


        

        REF_EPOC = char(epoc(streamAorC));
        ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
        ARTIFACT465 = Inf;% variable created for artifact removal for 465 store
        if streamAorC == 1
            if ~isfield(data.epocs, 'aRw_')
                disp('Epoc for position A is missing...')
                continue
            end
            disp("Plotting stream A")
            STREAM_STORE1 = 'x405A'; % name of the 405 store
            STREAM_STORE2 = 'x465A'; % name of the 465 store
            TITLE{batch} = strcat(animalIDA," ",taskA," ",EPOCNAME);
            savepath = strcat(figure_savepath,animalIDA,"/");
            if not(isfolder(savepath))
                mkdir(savepath);
            end
        elseif streamAorC == 2
            if ~isfield(data.epocs, 'bRw_')
                disp('Epoc for position C is missing...')
                continue
            end
            disp("Plotting stream C")
            STREAM_STORE1 = 'x405C'; % name of the 405 store
            STREAM_STORE2 = 'x465C'; % name of the 465 store
            TITLE{batch} = strcat(animalIDC," ",taskC," ",EPOCNAME);
            savepath = strcat(figure_savepath,animalIDC,"/");
            if not(isfolder(savepath))
                mkdir(savepath);
            end
            
        end
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
        
        % downsample and average 405 signal
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
        
        % Plot the 405 and 465 average signals
        figure(1);
        subplot(3,1,1);
        plot(ts1, meanSignal1, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 3); hold on;
        plot(ts2, meanSignal2, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); 
        
        % Plot vertical line at epoch onset, time = 0
        line([0 0], [(min(F465(:) - dcSignal2)), ((max(F465(:)) - dcSignal2))], 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)
        
        %Create the standard error bands for the 405 signal
        XX = [ts1, fliplr(ts1)];
        YY = [meanSignal1 + stdSignal1, fliplr(meanSignal1 - stdSignal1)];
        % 
        % Plot filled standard error bands.
        h = fill(XX, YY, 'g');
        set(h, 'facealpha',.25,'edgecolor','none')
        
        % Repeat for 465
        XX = [ts2, fliplr(ts2)];
        YY = [meanSignal2 + stdSignal2, fliplr(meanSignal2 - stdSignal2)];
        h = fill(XX, YY, 'r');
        set(h, 'facealpha',.25,'edgecolor','none')
        % 
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('mV', 'FontSize', 12)
        title(sprintf(TITLE{batch}, numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
        set(gcf, 'Position',[100, 100, 800, 500])
        
        % Heat Map based on z score of 405 fit subtracted 465
        % Fitting 405 channel onto 465 channel to detrend signal bleaching
        % Scale and fit data
        % Algorithm sourced from Tom Davidson's Github:
        % https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
        
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
        
        % Plot heat map
        subplot(3,1,2);
        imagesc(ts2, 1, zall);
        colormap('jet'); % c1 = colorbar; 
        title(sprintf('Z-Score Heat Map', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14);
        ylabel('Trials', 'FontSize', 12);
        
        %Fill band values for second subplot. Doing here to scale onset bar
        %correctly
        XX = [ts2, fliplr(ts2)];
        YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];
        
        subplot(3,1,3)
        plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
        line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
        
        h = fill(XX, YY, 'r');
        set(h, 'facealpha',.25,'edgecolor','none')
        
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('Z-score', 'FontSize', 12)
        title(sprintf('465 nm Z-Score', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
        c2 = colorbar;
        %%
        figure(2);
        plot(ts2, zall)
        line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
        
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('Z-score', 'FontSize', 12)
        title(sprintf('465 nm Z-Score', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
        
        %%
        fig{batch} = figure(3);
        subplot(2,3,[1,2,4,5])
        plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
        line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)
        if addVertLine == 1
            line([vertLineX vertLineX], [min(YY*1.5), max(YY*1.5)], 'Color', [1 0 0], 'LineWidth', 2)
        end
        
        h = fill(XX, YY, 'b');
        set(h, 'facealpha',.25,'edgecolor','none')
        
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',18)
        ylabel('Z-score +/- SEM', 'FontSize', 18)
        title(sprintf(TITLE{batch}, numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 18)
        box off
        
        subplot(2,3,6);
        imagesc(ts2, 1, zall);
        colormap('jet'); colorbar; 
        title(sprintf('Z-Score/Trial', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 16);
        xlabel('Time, s', 'FontSize', 12);
        ylabel('Trial', 'FontSize', 12);
        
        % Fill band values for second subplot. Doing here to scale onset bar
        % correctly
        XX = [ts2, fliplr(ts2)];
        YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];
        
        savdisp = strcat("Saving figure for: ",TITLE{batch});
        disp(savdisp)
        
        filename = strcat(savepath,TITLE{batch},savetype);
        print(filename,'-dpdf','-bestfit');
        
        
        clf reset   
    end
end
close all
if length(myFiles) > 1
    fprintf("Finished plotting and saving %d figures...\n",length(myFiles))
elseif length(myFiles) == 1
    fprintf("Finished plotting and saving %d figure...\n",length(myFiles))
end

NERD_STATS(toc,length(myFiles));

