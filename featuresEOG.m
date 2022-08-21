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
% correlates with values near -1.0, whereas NREM sleep correlates
% with values near +1.0.
% strongly correlates with values near -1.0. NREM sleep does not 
% exhibit this kind of behavior.
% The second column contains the sleep stage Annotation of that 
% window.
% -------------------------------------------------------------------

function [features] = featuresEOG(X)
    N = size(X,1);

    % Initialize an empty table to store
    % the features and sleep stage Annotations
    types = ["double" "double" "double" "string"];
    names = ["low" "mid" "high" "Annotations"];
    
    features = table(               ...
        'Size',             [N 4],  ...
        'VariableTypes',    types,  ...
        'VariableNames',    names);

    features.Annotations = X.Annotations;

    % MRA decomposition on 1st EOG channel
    dwt1 = mraEEG(X,"EOGE1_M2");

    % MRA decomposition on 2nd EOG channel
    dwt2 = mraEEG(X,"EOGE2_M2");

    for i = 1:1:N
        % normalized DWT coefficients of 0Hz-4Hz scale
        low1 = cell2mat(dwt1{i,"delta"});
        low2 = cell2mat(dwt2{i,"delta"});
        low1 = (low1 - mean(low1)) / std(low1);
        low2 = (low2 - mean(low2)) / std(low2);

        % normalized DWT coefficients of 4Hz-8Hz scale
        mid1 = cell2mat(dwt1{i,"theta"});
        mid2 = cell2mat(dwt2{i,"theta"});
        mid1 = (mid1 - mean(mid1)) / std(mid1);
        mid2 = (mid2 - mean(mid2)) / std(mid2);

        % normalized DWT coefficients of 8Hz-16Hz scale
        high1 = cell2mat(dwt1{i,"alpha"});
        high2 = cell2mat(dwt2{i,"alpha"});
        high1 = (high1 - mean(high1)) / std(high1);
        high2 = (high2 - mean(high2)) / std(high2);

        % cross-correlation of EOG channels in different scales
        features{i,"low"}  = mean(low1 .* low2);
        features{i,"mid"}  = mean(mid1 .* mid2);
        features{i,"high"} = mean(high1 .* high2);
    end
end
