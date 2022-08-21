function [features] = featuresEOG(X)
    N = size(X,1);

    % Initialize an empty table to store
    % the features and sleep stage Annotations
    types = ["double" "string"];
    names = ["corr" "Annotations"];
    
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
        features{i,"corr"} = mean(x.*y);
    end
end