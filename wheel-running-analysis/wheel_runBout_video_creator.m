%% MATLAB Script: Overlaying a Scrolling Signal Under a Video
% This script reads an input .avi video and a 1D signal, then outputs a new
% video that shows the video frames in the top half and a scrolling signal plot 
% in the bottom half. In the signal plot, only the most recent 'windowLength' seconds
% of data are shown, with a vertical red line indicating the current time.
%
% USER PARAMETERS (adjust these as needed):
channel = 1;
N = 100; % downsample
videoFile = '/Volumes/CUDADRIVE/eCB Wheel Harrison & Alex/Videos/output_1463 Box 1_3-7-25.avi';         % Input video filename (.avi)
outputVideoFile = '/Volumes/CUDADRIVE/eCB Wheel Harrison & Alex/Videos/1463_1.avi';  % Output video filename
BLOCKPATH = '/Volumes/CUDADRIVE/eCB Wheel Harrison & Alex/tanks/1463_Wheel1_EMPTY Scored';

signalSampleRate = 1017/N;         % Signal sampling rate in Hz (adjust as needed)
windowLength = 10;                % Time window in seconds to display for the scrolling signal

data = TDTbin2mat(BLOCKPATH);
%Stream Stores%
if channel == 1
    ISOS = 'x405A'; % name of the 405A store
    SIGNAL = 'x465A'; % name of the 465A store
elseif channel == 2
    ISOS = 'x405C';
    SIGNAL = 'x465C';
else
    error('Unknown channel number')
end
%time array used for all streams%
time = (1:length(data.streams.(SIGNAL).data))/data.streams.(SIGNAL).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%
t = 20; % time threshold below which we will discard
ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(SIGNAL).data = data.streams.(SIGNAL).data(ind:end);
data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);
min_time = min(time);

%downsample streams and time array by N times%
data.streams.(ISOS).data = downsample(data.streams.(ISOS).data, N);
data.streams.(SIGNAL).data = downsample(data.streams.(SIGNAL).data, N);
time = downsample(time, N);

%detrend & dFF%
bls = polyfit(data.streams.(ISOS).data,data.streams.(SIGNAL).data,1);
Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
Y_dF_all = data.streams.(SIGNAL).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
signal = detrend(dFF);
signal = smoothdata(signal,'movmean',smoothFactor);

%% Set Up
% Create a time vector for the signal.
tSignal = (0:length(signal)-1) / signalSampleRate;

% Open the input video.
vidObj = VideoReader(videoFile);
frameRate = vidObj.FrameRate;

% Prepare the output video writer.
outputVideo = VideoWriter(outputVideoFile);
outputVideo.FrameRate = frameRate;
open(outputVideo);

% Create a figure for displaying the combined view.
% The figure is set to a fixed size. Adjust the 'Position' as needed.
figureHandle = figure('Visible','off','Position',[100, 100, 1000, 800]);

% Create two subplots:
%   - The top subplot will display the video frame.
%   - The bottom subplot will display the scrolling signal.
videoAx = subplot(2,1,1);
signalAx = subplot(2,1,2);

%% Process Video Frames and Update Plot
frameCount = 0;
while hasFrame(vidObj)
    frameCount = frameCount + 1;
    currentFrame = readFrame(vidObj);
    
    % Calculate the current video time (in seconds).
    currentTime = (frameCount - 1) / frameRate;
    
    % ----- Top Axes: Display Video Frame -----
    axes(videoAx);         % Switch to video axes
    imshow(currentFrame);  % Show the video frame
    title(sprintf('Video Frame %d - Time: %.2f sec', frameCount, currentTime));
    
    % ----- Bottom Axes: Plot the Scrolling Signal -----
    % Define the time window for the signal display.
    % Here we show the previous windowLength seconds up to the current time.
    tMin = currentTime - windowLength;
    tMax = currentTime;
    if tMax <= tMin
        tMax = tMin + eps;
    end
    
    % Find the indices of the signal that fall within the current time window.
    indices = find(tSignal >= tMin & tSignal <= tMax);
    
    % Plot the portion of the signal.
    axes(signalAx);    % Switch to signal axes
    plot(tSignal(indices), signal(indices), 'b-', 'LineWidth', 1.5);
    hold on;
    % Draw a vertical red dashed line at the current time to indicate the "live" moment.
    yLimits = ylim;
    plot([currentTime currentTime], yLimits, 'r--', 'LineWidth', 2);
    hold off;
    xlabel('Time (s)');
    ylabel('Signal');
    title('Scrolling Signal');
    xlim([tMin, tMax]);
    
    % Force the figure to update (without displaying it).
    drawnow;
    
    % Capture the current figure as a frame.
    frameData = getframe(figureHandle);
    
    % Write the captured frame to the output video.
    writeVideo(outputVideo, frameData);
end

%% Finalize
close(outputVideo);       % Close the video file.
close(figureHandle);        % Close the figure.

disp('Video creation complete.');