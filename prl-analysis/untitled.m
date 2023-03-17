% Generate some example data
signal_data = randn(10, 100); % 10 rows of 100 samples each
associated_data = randn(10, 2); % 10 rows of 2 associated values each

% Create a new figure
figure;
hold on;

% Loop through each row of signal data and plot the associated data as arrows
for i = 1:size(signal_data, 1)
    % Extract the ith row of signal data and associated data
    signal_row = signal_data(i, :);
    associated_row = associated_data(i, :);
    
    % Plot the signal data
    plot(signal_row);
    
    % Plot the first arrow in red
    quiver(1, associated_row(1), signal_row(1) - associated_row(1), 0, 'r');
    
    % Loop through the rest of the associated data and plot arrows in green
    for j = 2:length(associated_row)
        x_pos = find(signal_row == associated_row(j));
        if ~isempty(x_pos)
            quiver(x_pos, associated_row(j), signal_row(x_pos) - associated_row(j), zeros(size(x_pos)), 'g');
        end
    end
end

% Add axis labels and title
xlabel('Sample Number');
ylabel('Associated Value');
title('Signal Data with Associated Data');

% Hold off to end plotting for current figure
hold off;


