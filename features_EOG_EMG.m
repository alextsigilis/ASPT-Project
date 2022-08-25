% ===================================================================
%
% Author: Chrysa Doulou, Christodoulos Michaelides
% Date: August 21st, 2022
% -------------------------------------------------------------------
%
% Function Description: 
% EOG signals provide useful information for distinguishing between 
% REM and Non-REM sleep. Chin-EMG recordings are also useful for 
% distinguishing between the first stages of sleep and wakefulness.
% This function extracts features from such recordings which can then
% be used by an automatic sleep-stage scoring system.
% -------------------------------------------------------------------
%
% Arguments List: (X)
%
% X: (table) A table which contains EEG,EOG,EMG and ECG recordings.
% You should use the loadEDF function to obtain this table.
% -------------------------------------------------------------------
%
% Return List: (features)
%
% features: (table) a table with two columns. 
% The 1st column contains the cross-correlation coefficient 
% of the EOG signals estimated in a 30sec window. REM sleep 
% correlates with values near -1.0, whereas NREM sleep strongly
% correlates with values near +1.0. In other words, the
% cross-correlation is a very useful tool for detecting conjugate
% waveforms in EOG signals, which are a common occurence
% during REM sleep.
% The 2nd and 3rd columns contain the total energy of the EOG channels 
% inside the 0.10Hz-1.00Hz frequency band (ECB).
% The 4th column contains the standard deviation of the chin-EMG.
% Small values indicate muscle relaxiation and loss of consciousness,
% whereas large values indicate wakefulness and alertness.
% The 5th column contains the peak-to-peak voltage of the EMG 
% recordings.
% The 6th column contains the sleep stage Annotations.
%
% useDWT: (boolean) by setting useDWT to true, the discrete wavelet
% transform is used to decompose the EOG signals into five 
% frequency scales. The cross-correlation of the EOG channels is
% estimated by using the DWT coefficients of the 0Hz-4Hz frequency 
% scale. This method can lead to better results because it 
% removes high-frequency noise.
% -------------------------------------------------------------------

function [features] = features_EOG_EMG(X, useDWT)
    % Number of 30sec epochs
    N = size(X,1);
    if N == 0 error("X is empty"); end

    % Initialize an empty table to store
    % features and sleep stage Annotations
    types = [               ...
        "double",           ...
        "double",           ...
        "double",           ...
        "double",           ...
        "double",           ...
        "string"];
    
    names = [               ...
        "xcorr",            ...         % EOG cross-correlation
        "ECB1",             ...         % ECB (1st EOG channel)
        "ECB2",             ...         % ECB (2nd EOG channel)
        "stdEMG",           ...         % chin-EMG standard deviation
        "VppEMG",           ...         % chin-EMG peak to peak voltage
        "Annotations"];                 % sleep stage Annotations
    
    features = table(               ...
        'Size',             [N 6],  ...
        'VariableTypes',    types,  ...
        'VariableNames',    names);

    % Copy sleep stage Annotations 
    features.Annotations = X.Annotations;

    % T:   (float) duration of EOG recordings in seconds
    % fs:  (float) sampling frequency of EOG recordings
    % M:   (integer) number of samples in EOG recordings
    % f:   (1D array) frequency axis of EOG power spectrum
    % idx: (1D array) vector of indices for frequencies 
    %       between 0.1Hz and 1.0Hz from the EOG power spectra.
    
    % idx = 4:1:31;
    idx = 11:1:127;
    % ------------------------------------------------------- 
    % T = 30; fs = 256; M = fs * T;
    % f = (-fs/2):(fs/M):(+fs/2-fs/M); f = ifftshift(f);
    % idx = find(f >= 0.3 & f <= 0.5);
    % -------------------------------------------------------

    for i = 1:1:N
        % Extract 30sec segments from
        % the EOG and EMG recordings.
        x = cell2mat(X{i,"EOGE1_M2"});
        y = cell2mat(X{i,"EOGE2_M2"});
        z = cell2mat(X{i,"EMGChin"});
        
        % ECB values for EOG channels
        Xf = fft(x) / numel(x); Xf = abs(Xf(idx)).^2; 
        Yf = fft(y) / numel(y); Yf = abs(Yf(idx)).^2;
        features{i,"ECB1"} = sum(Xf) / numel(Xf);
        features{i,"ECB2"} = sum(Yf) / numel(Yf);
        
        % Wavelet Denoising
        if useDWT == true
            % wavelet decomposition
            [x, lx] = wavedec(x, 5, 'db5');
            [y, ly] = wavedec(y, 5, 'db5');

            % DWT approximation
            % coefficients (0Hz-4Hz)
            x = appcoef(x, lx, 'db5', 5);
            y = appcoef(y, ly, 'db5', 5);
        end

        % EMG standard deviation
        features{i,"stdEMG"} = std(z);

        % EMG peak-to-peak voltage
        features{i,"VppEMG"} = max(z) - min(z);

        % EOG cross-correlation
        epsilon = 1e-4;              % numerical stability
        mx = mean(x); sx = std(x);
        my = mean(y); sy = std(y);
        mxy = mean(x.*y); 
        features{i,"xcorr"} = (mxy-mx*my)/(sx*sy+epsilon);
    end
end
