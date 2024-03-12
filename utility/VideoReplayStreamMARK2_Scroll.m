%% Video Replay with Streaming Data
%
%  Import video and stream data into Matlab
%  Plot the video and stream
%  Good for replay analysis

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clear all; clc;
% SDKPATH = 'C:\Users\Andrew Villa\Desktop\Matlab\WheelFP\PhotometryScripts\TDTSDK\TDTbin2mat'; % or whatever path you extracted the SDK zip into
% addpath(genpath(SDKPATH));
BLOCKPATH = '\Users\Andrew Villa\Desktop\Matlab\WheelFP\Wheel-DA13_3';
data = TDTbin2mat(BLOCKPATH, 'TYPE',{'epocs', 'scalars', 'streams'}); %,'T1',1500,'T2',1600,
savePATH = '/Users/Andrew Villa/Desktop/Matlab';
%% Importing the Data
% This example assumes you downloaded our
% <https://www.tdt.com/files/examples/TDTExampleData.zip example data sets>
% and extracted it into the \TDTMatlabSDK\Examples\ directory. To import your own data, replace
% 'BLOCKPATH' with the path to your own data block.
%
% In Synapse, you can find the block path in the database. Go to Menu --> History. 
% Find your block, then Right-Click --> Copy path to clipboard.
% BLOCKPATH           = fullfile(DATAPATH,'Subject1-211115-094936');

%Stream Stores% Set stream store variables to create objects
DLS_ISOS = 'x405A';     % name of the 405A store DLS
DLS_DA = 'x465A';       % name of the 465A store DLS
NAc_ISOS = 'x405C';     % name of the 405C store NAc
NAc_DA = 'x465C';       % name of the 465C store NAc
% EPOC = 'runBout';       % name of running bout epoc created later%
VID_STORE = 'Cam2';     % video store name
START_FRAME = 29980;    % first frame index (use 1 for beginning of video)
END_FRAME = 31980;      % last frame index (use -1 for end of video)
ROLLING = 10;          % rolling window, in seconds (use -1 for none)
CREATE_OUTPUT_VIDEO = 1;% set to 0 to skip writing the output video
VIDEO_OUTPUT_PATH = BLOCKPATH; % where the output video should go
data = TDTbin2mat(BLOCKPATH); %set TDT bin2mat path

%% %% %% epocs created in TDT OpenScope save in the notes of Cam1 %% 

% %combine index with timestamp data from Cam1 notes%
ind = double(data.epocs.Cam1.notes.index);
ts = data.epocs.Cam1.notes.ts;
var1 = [ind ts];


%average around every Nth point and downsample Nx% Downsample data
N = 10; % multiplicative for downsampling
data.streams.(NAc_DA).data = arrayfun(@(i)...
    mean(data.streams.(NAc_DA).data(i:i+N-1)),...
    1:N:length(data.streams.(NAc_DA).data)-N+1);
data.streams.(NAc_ISOS).data = arrayfun(@(i)...
    mean(data.streams.(NAc_ISOS).data(i:i+N-1)),...
    1:N:length(data.streams.(NAc_ISOS).data)-N+1);
data.streams.(DLS_DA).data = arrayfun(@(i)...
    mean(data.streams.(DLS_DA).data(i:i+N-1)),...
    1:N:length(data.streams.(DLS_DA).data)-N+1);
data.streams.(DLS_ISOS).data = arrayfun(@(i)...
    mean(data.streams.(DLS_ISOS).data(i:i+N-1)),...
    1:N:length(data.streams.(DLS_ISOS).data)-N+1);

%establish sample rate variables %decimate time array and match length to demodulated stream%
fs = data.streams.(DLS_DA).fs/N; %DLS% /N as recall we downsampled by N = 10 earlier
fs2 = data.streams.(NAc_DA).fs/N;%NAc%

%Calculations for data to plot%
% Create mean signal, standard error of signal, and DC offset of DLS_ISOS signal
meanSignal1 = mean(DLS_ISOS);
stdSignal1 = std(double(DLS_ISOS))/sqrt(size(DLS_ISOS,1));
dcSignal1 = mean(meanSignal1);
medianSignal1 = median(DLS_ISOS); % median of fluorescence signal for dF/F calc

% Create mean signal, standard error of signal, and DC offset of DLS_DA signal
meanSignal1 = mean(DLS_DA);
stdSignal1 = std(double(DLS_DA))/sqrt(size(DLS_DA,1));
dcSignal1 = mean(meanSignal1);
medianSignal1 = median(DLS_DA); % median of fluorescence signal for dF/F calc

% Create mean signal, standard error of signal, and DC offset of NAc_ISOS signal
meanSignal1 = mean(NAc_ISOS);
stdSignal1 = std(double(NAc_ISOS))/sqrt(size(NAc_ISOS,1));
dcSignal1 = mean(meanSignal1);
medianSignal1 = median(NAc_ISOS); % median of fluorescence signal for dF/F calc

% Create mean signal, standard error of signal, and DC offset of NAc_DA signal
meanSignal1 = mean(NAc_DA);
stdSignal1 = std(double(NAc_DA))/sqrt(size(NAc_DA,1));
dcSignal1 = mean(meanSignal1);
medianSignal1 = median(NAc_DA); % median of fluorescence signal for dF/F calc

%set index variables for Zall to match video data
ct = 1; 
k = START_FRAME:END_FRAME;
recording_ts = data.epocs.(VID_STORE).onset(k);
 stream_ind = round(recording_ts * data.streams.(DLS_DA).fs);
    if ct == 1
        start_ind = stream_ind;
    end

%detrend & dFF & Z-All% Calculations for creating Zall values to plot
%DLS%
%465/405A


bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF

zall = zeros(size(Y_dF_all));
for i = 1:size(Y_dF_all,1)
    ind = ts;%(1,:);% index %some errors may occur here due to index array restrictions
    zb = mean(Y_dF_all); % baseline mean
    zsd = std(Y_dF_all); % baseline stdev
    zall=(Y_dF_all - zb)/zsd; % Z score per bin
end
zall = detrend(zall); %detrend zall

%NAc%
%465/405C%
bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF


zall2 = zeros(size(Y_dF_all2)); %NAc
for i = 1:size(Y_dF_all2,1)
    ind = ts;%(1,:);    %index % some errors may occur here due to index array restrictions
    zb2 = mean(Y_dF_all2); % baseline period mean (-10sec to -6sec)
    zsd2 = std(Y_dF_all2); % baseline period stdev
    zall2=(Y_dF_all2 - zb2)/zsd2; % Z score per bin
end
zall2 = detrend(zall2); %detrend zall

%%
% Read video file.
vvv = dir([BLOCKPATH filesep '*' VID_STORE '.avi']);
vid_filename = [vvv.folder filesep vvv.name];
fprintf('reading file %s\n', vid_filename);
myvideo = VideoReader(vid_filename);

%%
% Get data specs.
max_frames = length(data.epocs.(VID_STORE).onset);
if END_FRAME < 1
    END_FRAME = max_frames;
end
max_ts = data.epocs.(VID_STORE).onset(end);
expected_fps = max_frames / max_ts;
max_x = max(size(data.streams.(DLS_DA).data));

%%
% Make array of images if we're outputting a video.
if CREATE_OUTPUT_VIDEO
     M(END_FRAME-START_FRAME+1) = struct('cdata',[],'colormap',[]);
end

%%
% Create figure.
h = figure;
h.Position = [500 500 560 560];

%%
% The main loop.

tic
ct = 1;
for k = START_FRAME:END_FRAME
    % grab one image
    im = read(myvideo, k);
    
    subplot(4,1,[1 2])

    % plot it
    image(im)
    if ct == 1
        % hide x and y pixel axes
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'nextplot','replacechildren')
    end
    recording_ts = data.epocs.(VID_STORE).onset(k);

    % set title
    title_text = sprintf('%s frame %d of %d, t = %.2fs', VID_STORE, k, END_FRAME, recording_ts);
    title(title_text);
    
    % plot stream DLS in another subplot
    subplot(4,1,3)

    stream_ind = (round(recording_ts * data.streams.(DLS_DA).fs)/N);
    if ct == 1
        start_ind = stream_ind;
        end_ind = (round(data.epocs.(VID_STORE).onset(END_FRAME) * data.streams.(DLS_DA).fs)/N);
    end
    
    if ROLLING > 0
%         stream_ind = round(recording_ts * data.streams.(DLS_DA).fs);
        max_ts = data.epocs.(VID_STORE).onset(k);
        start_ind = (max(round(max(recording_ts - ROLLING, 0) * data.streams.(DLS_DA).fs), 1)/N);
        t = max_ts .* (start_ind:stream_ind) / stream_ind;
    else
        t = max_ts .* (start_ind:stream_ind) ./ max_x;
    end
    
    plot(t, zall(start_ind:stream_ind), 'b', 'LineWidth', 2) %data.streams.(DLS_DA).data
    if ct == 1
        t1 = max_ts .* start_ind / max_x;
        t2 = max_ts .* end_ind / max_x;
        y1 = 0;
        y2 = max(zall(start_ind:end_ind)); %data.streams.(DLS_DA).data
        axis([t1, t2, y1, y2])
        grid on;
        title('DLS DA', 'color', 'b')
        xlabel('time, s')
        ylabel('Z-score')
        ylim manual
        set(gca,'nextplot','replacechildren') % maintains the axis properties next time, improves speed
    end
    if ROLLING > 0
        t1 = t(1);
        t2 = max(t(end), ROLLING);
        axis([t1, t2, y1, y2]);
    end
    axis tight;
    % force the plot to update
    drawnow;
    
%     if CREATE_OUTPUT_VIDEO
%         M(ct) = getframe(gcf); % get the whole figure
%     end
% 
%     % slow down to match video fps
%     expected_el = ct / expected_fps;
%     ddd = expected_el - toc;
%     if ddd > 0, pause(ddd); end
%     ct = ct + 1;
% end

    % plot stream NAc in another subplot
    subplot(4,1,4)

    stream_ind = (round(recording_ts * data.streams.(NAc_DA).fs)/N);
    if ct == 1
        start_ind = stream_ind;
        end_ind = (round(data.epocs.(VID_STORE).onset(END_FRAME) * data.streams.(NAc_DA).fs)/N);
    end
    
    if ROLLING > 0
%         stream_ind = round(recording_ts * data.streams.(NAc_DA).fs);
        max_ts = data.epocs.(VID_STORE).onset(k);
        start_ind = (max(round(max(recording_ts - ROLLING, 0) * data.streams.(NAc_DA).fs), 1)/N);
        t = max_ts .* (start_ind:stream_ind) / stream_ind;
    else
        t = max_ts .* (start_ind:stream_ind) ./ max_x;
    end
    
    plot(t, zall2(start_ind:stream_ind), 'r', 'LineWidth', 2) %data.streams.(NAc_DA).data
    if ct == 1
        t1 = max_ts .* start_ind / max_x;
        t2 = max_ts .* end_ind / max_x;
        y1 = 0;
        y2 = max(zall2(start_ind:end_ind)); %data.streams.(NAc_DA).data
        axis([t1, t2, y1, y2])
        grid on;
        title('NAc DA', 'color', 'r')
        xlabel('time, s')
        ylabel('Z-Score')
        ylim manual
        set(gca,'nextplot','replacechildren') % maintains the axis properties next time, improves speed
    end
    if ROLLING > 0
        t1 = t(1);
        t2 = max(t(end), ROLLING);
        axis([t1, t2, y1, y2]);
    end
axis tight;
    % force the plot to update
    drawnow;
    
    if CREATE_OUTPUT_VIDEO
        M(ct) = getframe(gcf); % get the whole figure
    end

    % slow down to match video fps
    expected_el = ct / expected_fps;
    ddd = expected_el - toc;
    if ddd > 0, pause(ddd); end
    ct = ct + 1;
end

disp('done playing')

%%
% Create the output video file of figure with same FPS as original.
if CREATE_OUTPUT_VIDEO
    out_file = [VIDEO_OUTPUT_PATH filesep strrep(vvv.name, '.avi', '_output.avi')];
    fprintf('writing video file %s\n', out_file);
    out_video = VideoWriter(out_file);
    out_video.FrameRate = expected_fps;
    open(out_video);
    for k = 1:length(M)
        writeVideo(out_video, M(k));
    end
    close(out_video)
end
