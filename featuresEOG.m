% ===================================================================
%
% Author: Christodoulos Michaelides
% Date: August 21st, 2022
% -------------------------------------------------------------------
%
% Function Description: 
% Distinguishing between REM and NREM sleep usually requires
% visual inspection of two EOG recordings. This functions extracts
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
% The first column contains the cross-correlation coefficient 
% of the EOG signals estimated in a 30sec window. REM sleep 
<<<<<<< HEAD
% correlates with values near -1.0, whereas NREM sleep correlates
% with values near +1.0.
=======
% strongly correlates with values near -1.0. NREM sleep does not 
% exhibit this kind of behavior.
>>>>>>> 74d490bc7bd9fa82c6cc548781c86adc8e213a6f
% The second column contains the sleep stage Annotation of that 
% window.
% -------------------------------------------------------------------

function [features] = featuresEOG(X)
    N = size(X,1);

    % Initialize an empty table to store
    % the features and sleep stage Annotations
    types = ["double" "string"];
    names = ["xcorr" "Annotations"];
    
    features = table(               ...
        'Size',             [N 2],  ...
        'VariableTypes',    types,  ...
        'VariableNames',    names);

    % Copy the sleep stage Annotations
    features.Annotations = X.Annotations;

    for i = 1:1:N
        % Extract two 30sec epochs 
        % from the EOG channels
        x = cell2mat(X{i,"EOGE1_M2"});
        y = cell2mat(X{i,"EOGE2_M2"});
        
        % Normalize the samples
        x = (x - mean(x)) / std(x);
        y = (y - mean(y)) / std(y);
        
        % Estimate the cross-correlation coefficient
        features{i,"xcorr"} = mean(x .* y);
    end
end
