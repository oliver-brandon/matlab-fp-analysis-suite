function make_toy_photometry_session(outFile)
% make_toy_photometry_session Create a small representative TDTbin2mat-like struct.
% Usage:
%   make_toy_photometry_session('toy_photometry_session.mat')

if nargin < 1
    outFile = 'toy_photometry_session.mat';
end

rng(0);

% --- Sampling and duration ---
fs = 1017;              % Hz
dur_s = 30;             % seconds (small)
n = dur_s * fs;
t = (0:n-1) / fs;

% --- Create a shared "motion" artifact component ---
motion = 0.02 * sin(2*pi*1.3*t) + 0.01 * randn(1, n);

% --- Iso (405) signals per site ---
x405A = 1.0 + 0.02*sin(2*pi*0.2*t) + motion + 0.005*randn(1, n);
x405C = 1.0 + 0.02*sin(2*pi*0.23*t + 0.4) + motion + 0.005*randn(1, n);

% --- "True" neural events (for 465/560) ---
event_onsets = [5, 10, 15, 20, 25]; % seconds
kernel_t = (0:1/fs:2.0);
kernel = exp(-kernel_t/0.4) .* (1 - exp(-kernel_t/0.05)); % quick rise, slow decay
kernel = kernel / max(kernel);

ev = zeros(1, n);
for k = 1:numel(event_onsets)
    idx0 = round(event_onsets(k) * fs) + 1;
    idx1 = min(n, idx0 + numel(kernel) - 1);
    ev(idx0:idx1) = ev(idx0:idx1) + 0.05 * kernel(1:(idx1-idx0+1));
end

% --- Signal channels: correlated with iso + added event component ---
% These are intentionally constructed so that regression/correction has something to do.
x465A = 0.8*x405A + 0.2 + ev + 0.006*randn(1, n);
x560A = 0.75*x405A + 0.25 + 0.6*ev + 0.006*randn(1, n);

x465C = 0.82*x405C + 0.18 + 0.8*ev + 0.006*randn(1, n);
x560C = 0.78*x405C + 0.22 + 0.5*ev + 0.006*randn(1, n);

% --- Build `data` struct in the style you rely on ---
data = struct();

% Streams: each store has `.data` and `.fs`
data.streams = struct();
data.streams.x405A = struct('data', x405A', 'fs', fs);
data.streams.x465A = struct('data', x465A', 'fs', fs);
data.streams.x560A = struct('data', x560A', 'fs', fs);

data.streams.x405C = struct('data', x405C', 'fs', fs);
data.streams.x465C = struct('data', x465C', 'fs', fs);
data.streams.x560C = struct('data', x560C', 'fs', fs);

% Epocs: at minimum `.onset` (seconds)
data.epocs = struct();
data.epocs.Stim = struct('onset', event_onsets(:)); % example ref epoc name "Stim"

% Info: arbitrary metadata
data.info = struct();
data.info.subject = 'toy_mouse_01';
data.info.session = 'toy_session';
data.info.fs_note = 'All streams fs=1017';
data.info.created = datestr(now);

save(outFile, 'data', '-v7'); % classic .mat for maximal compatibility
fprintf('Wrote %s\n', outFile);
end
