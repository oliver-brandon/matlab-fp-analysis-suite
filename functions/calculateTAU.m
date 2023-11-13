function tau = calculateTAU(signal, time)
    % This function calculates the decay constant tau for a fiber photometry
    % signal.
    % 
    % Inputs:
    % - time: Vector containing time points corresponding to the signal
    % - signal: Vector containing the fiber photometry signal
    %
    % Output:
    % - tau: Calculated decay constant (tau)
    
    signal = signal';
    time = time';
    % Fit options
    fit_options = fitoptions('Method', 'NonlinearLeastSquares', ...
                             'Lower', [0, 0, 0], ...
                             'Upper', [Inf, Inf, Inf], ...
                             'StartPoint', [1, 1, 1]);
    
    % Define fit type (exponential decay)
    fit_type = fittype('A*exp(-x/tau) + C', ...
                       'options', fit_options);
    
    % Perform the fit
    [fitted_model, ~] = fit(time(:), signal(:), fit_type);
    
    % Extract tau from the fitted model parameters
    tau = fitted_model.tau;
end
