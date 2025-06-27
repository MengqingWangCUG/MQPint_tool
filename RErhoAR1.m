function [rho] = RErhoAR1(data)
    % RHOAR1 calculates the lag-1 autocorrelation coefficient for the given time series.
    % This function assumes that the data is resampled with a constant sampling step.
    %
    % Inputs:
    %   data - A column vector containing the time series data (with constant sampling step).
    %
    % Outputs:
    %   rho - The lag-1 autocorrelation coefficient of the time series.
    %
    % Author: Doroth¨¦e Husson, November 2012
    % Modified for code clarity and optimization

    % Length of the input data
    n = length(data);
    
    % Check if data is valid
    if n < 2
        error('Data must contain at least two values.');
    end
    
    % Compute the mean of the data
    meanData = mean(data);
    
    % Detrend the data by subtracting the mean
    detrendedData = data - meanData;
    
    % Compute the lag-1 autocorrelation coefficient (rho)
    % Sum of the products of consecutive data points
    numerator = sum(detrendedData(2:end) .* detrendedData(1:end-1));
    
    % Sum of squared values of the previous data points
    denominator = sum(detrendedData(1:end-1).^2);
    
    % Compute the lag-1 autocorrelation coefficient
    rho = numerator / denominator;
end
