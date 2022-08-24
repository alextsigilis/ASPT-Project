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
% The 2nd and 3rd columns contain the total energy of the EOG
% channels inside the 0.35Hz-0.50Hz frequency band (ECB). Those 
% frequencies appear during rapid eye movements and they are a 
% a strong indicator of REM sleep.
% The 4th column contains the standard deviation of the chin-EMG.
% Small values indicate muscle relaxiation and loss of consciousness,
% whereas large values indicate wakefulness and alertness.
% The 5th column contains the sleep stage Annotations.
%
% useDWT: (boolean) by setting useDWT to true, the discrete wavelet
% transform is used to decompose the EOG and EMG signals into five 
% frequency scales. The standard deviation of the EMG and the
% cross-correlation of the EOG channels are estimated by using the 
% DWT coefficients of the 0Hz-4Hz frequency scale. This method 
% can lead to better results because it removes high-frequency noise. 
% The estimation of the ECB values is not affected by useDWT.
% -------------------------------------------------------------------

function [features] = features_EOG_EMG(X, useDWT)
    % Number of 30sec epochs
    N = size(X,1);

    % Initialize an empty table to store
    % the features and sleep stage Annotations
    types = [               ...
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
        "Annotations"];                 % sleep stage Annotations
    
    features = table(               ...
        'Size',             [N 5],  ...
        'VariableTypes',    types,  ...
        'VariableNames',    names);

    % Copy sleep stage Annotations 
    features.Annotations = X.Annotations;

    for i = 1:1:N
        % Extract 30sec segments from
        % the EOG and EMG recordings.
        x = X{i,"EOGE1_M2"};
        y = X{i,"EOGE2_M2"};
        z = X{i,"EMGChin"};
        
        % ECB values for EOG channels
        features{i,"ECB1"} = nan;
        features{i,"ECB2"} = nan;
        
        % Wavelet Denoising
        if useDWT == true
            % wavelet decomposition
            [x, lx] = wavedec(x, 5, 'db5');
            [y, ly] = wavedec(y, 5, 'db5');
            [z, lz] = wavedec(z, 5, 'db5');

            % DWT approximation
            % coefficients (0Hz-4Hz)
            x = appcoef(x, lx, 'db5', 5);
            y = appcoef(y, ly, 'db5', 5);
            z = appcoef(z, lz, 'db5', 5);
        end

        % EMG standard deviation
        features{i,"stdEMG"} = std(z);

        % EOG cross-correlation
        mx  = mean(x); my = mean(y);
        sx  = std(x); sy = std(y);
        mxy = mean(x.*y); 
        features{i,"xcorr"} = (mxy - mx*my)/(sx*sy);
    end
end
