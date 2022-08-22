% ===================================================================
%
% Author: Chrysa Doulou, Christodoulos Michaelides
% Date: August 21st, 2022
% -------------------------------------------------------------------
%
% Function Description: 
% Distinguishing between REM and NREM sleep usually requires
% visual inspection of EOG recordings. This function extracts
% useful features from the EOG channels which can then be used to  
% automate this decision process. 
% -------------------------------------------------------------------
%
% Arguments List: (X)
%
% X: (table) A table which contains EEG,EOG,EMG and ECG recordings.
% You should use loadEDF to obtain this table.
% -------------------------------------------------------------------
%
% Return List: (features)
%
% features: (table) a table with two columns. 
% The 1st column contains the cross-correlation coefficient 
% of the EOG signals estimated in a 30sec window. REM sleep 
% correlates with values near -1.0, whereas NREM sleep correlates
% with values near +1.0.
% The 2nd and 3rd column contain the standard deviation of the 
% EOG channels. Large variation corresponds to strong eye
% movements.
% The 4th and 5th columns contain the AUC values (area under
% curve) of the EOG channels. Large AUC values indicate strong
% eye movement.
% The 6th column contains the sleep stage Annotation of the 
% EOG signals.
% -------------------------------------------------------------------

function [features] = featuresEOG(X)
    % Number of 30sec epochs
    N = size(X,1);

    % Initialize an empty table to store
    % the features and sleep stage Annotations
    types = [               ...
        "double"            ...
        "double"            ...
        "double"            ...
        "double"            ...
        "double"            ...
        "string"];
    
    names = [               ...
        "xcorr"             ...
        "stdDev1"           ...
        "stdDev2"           ...
        "AUC1"              ...
        "AUC2"              ...
        "Annotations"];
    
    features = table(               ...
        'Size',             [N 6],  ...
        'VariableTypes',    types,  ...
        'VariableNames',    names);

    features.Annotations = X.Annotations;

    % Wavelet Decomposition
    dwt1 = mraEEG(X,"EOGE1_M2");
    dwt2 = mraEEG(X,"EOGE2_M2");

    for i = 1:1:N
        % Extract 30sec segments from the EOG signals
        x = cell2mat(dwt1{i,"delta"});
        y = cell2mat(dwt2{i,"delta"});

        % Calculate (normalized) AUC values
        features{i,"AUC1"} = sum(abs(x)) / numel(x);
        features{i,"AUC2"} = sum(abs(y)) / numel(y);

        % Calculate standard deviation
        features{i,"stdDev1"} = std(x);
        features{i,"stdDev2"} = std(y);

        % Calculate the cross correlation of the EOG channels
        w = (x - mean(x)) ./ std(x);
        z = (y - mean(y)) ./ std(y);
        features{i,"xcorr"} = mean(w .* z);
    end
end
