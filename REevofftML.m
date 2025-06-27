function [powRatio] = REevofftML(data, filename, fmin, fmax, window, fterm, nw)
    % This function calculates the evolutionary power spectral density using 
    % the multitaper method (MTM) and plots the ratio of power in the selected 
    % frequency band to the total power.
    %
    % Inputs:
    %   data    - A matrix where the first column is time and the second column 
    %             is the data values (evenly spaced time series).
    %   filename - A string representing the file name.
    %   fmin    - Minimum frequency for the power calculation (Hz).
    %   fmax    - Maximum frequency for the power calculation (Hz).
    %   window  - Window length for the sliding window analysis (s).
    %   fterm   - Upper frequency limit for the total power (Hz).
    %   nw      - The number of tapers to use in the multitaper method.
    %
    % Outputs:
    %   powRatio - The ratio of selected frequency power to total power.
    
    % Time step and Nyquist frequency calculation
    dt = mean(diff(data(:, 1)));  % Mean time step
    nyquist = 1 / (2 * dt);       % Nyquist frequency
    
    % Data size and preparation
    [nrow, ncol] = size(data);
    xdata = data(:, 2);           % Extract the data values (second column)
    npts = round(window / dt);    % Number of data points per window
    m = nrow - npts + 1;          % Number of windows for the analysis
    
    % Initialize matrices for power calculations
    pow = zeros(1, m);            % Power in the selected frequency band
    powAll = zeros(1, m);         % Total power in the frequency range [0, fterm]
    ratio = zeros(1, m);          % Power ratio (selected / total)
    
    % Results matrix for the multitaper method
    ss = zeros(m, length(pmtm(xdata, nw)));  % Store MTM results for each window
    
    % Perform the evolutionary PMTM calculations using a sliding window
    for i = 1:m
        windowData = xdata(i:i + npts - 1);  % Extract data for the current window
        windowData = detrend(windowData);    % Detrend the window data
        
        % Calculate the PMTM for the current window
        [p, ~] = pmtm(windowData, nw);
        ss(i, :) = p;  % Store the PMTM result
        
        % Frequency power calculation
        nfmin = max(1, ceil(length(p) * fmin / nyquist));  % Minimum frequency index
        nfmax = min(fterm, floor(length(p) * fmax / nyquist));  % Maximum frequency index
        totalPower = sum(p(nfmin:nfmax));  % Total power in the selected frequency range
        pow(i) = totalPower;
        
        % Total power up to fterm
        totalPowerAll = sum(p(1:fterm)); 
        powAll(i) = totalPowerAll;
        
        % Calculate power ratio
        ratio(i) = totalPower / totalPowerAll;
    end
    
    % Create time grid for plotting
    timeGrid = linspace(data(1, 1) + window / 2, data(end, 1) - window / 2, m);
    
    % Plot the results
    figure;
    
    subplot(3, 1, 1);
    plot(timeGrid, ratio);
    title(['Power Ratio (', num2str(fmin), ' Hz to ', num2str(fmax), ' Hz)']);
    xlabel('Time');
    ylabel('Power Ratio');
    
    subplot(3, 1, 2);
    plot(timeGrid, pow);
    title(['Total Power in ', num2str(fmin), ' Hz to ', num2str(fmax), ' Hz']);
    xlabel('Time');
    ylabel('Power');
    
    subplot(3, 1, 3);
    plot(timeGrid, powAll);
    title(['Total Power from 0 Hz to ', num2str(fterm), ' Hz']);
    xlabel('Time');
    ylabel('Total Power');
    
    % Set figure title
    suptitle(['Power Spectrum for ', filename]);
    
    % Return power ratio as output
    powRatio = ratio;
end
