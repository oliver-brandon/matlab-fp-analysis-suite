function numberedChoice = convertStringChoice(input)
    % Validate input
    if ~iscolumn(input) || ~iscellstr(input)
        error('Input must be a column of strings.');
    end
    
    % Preallocate output array
    numberedChoice = zeros(size(input));
    
    % Convert 'C' to 1 and 'I' to 2
    for i = 1:length(input)
        switch input{i}
            case 'C'
                numberedChoice(i) = 1;
            case 'I'
                numberedChoice(i) = 0;
            otherwise
                error('Input column contains invalid entries. Only ''C'' and ''I'' are allowed.');
        end
    end
    
    % Convert output to column vector
    numberedChoice = numberedChoice(:);
end