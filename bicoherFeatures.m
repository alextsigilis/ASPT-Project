% =======================================================
% Author Christodoulos Michaelides
% Date: August 25th, 2022
% -------------------------------------------------------
%
% Function Description: 
% This function can be used to extract features from the
% bicoherence estimations of EEG signals such as:
%   1) bicoherence entropy
%   2) bicoherence squared entropy
%   3) bicoherence cubed entropy
% -------------------------------------------------------
%
% Arguments List: (X, freq)
% X: (table) a table with two columns. The first should
% contain the bicoherence estimations and the
% second should contain the sleep stage Annotations.
% You should use bicEEG to obtain this table.
% -------------------------------------------------------
%
% Return List: (Y)
%
% TODO: add text
% =======================================================

function [Y] = bicoherFeatures(X)
    % Number of EEG segments (30sec epochs)
    N = size(X,1);

    % Initialize an empty table to store
    % features from the bicoherence matrices.
    types = ["double", "double", "double", "string"];
    names = ["ent1", "ent2", "ent3", "Annotations"];

    Y = table(                      ...
        'Size',          [N 4],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    for i = 1:1:N
        % Extract the entire bicoherence
        % matrix for a 30sec EEG segment
        bic = cell2mat(X{i,1});
        
        % Estimate bicoherence entropies 
        p = bic(:).^1; 
        Y{i,"ent1"} = -sum(p(:).*log2(p(:)),'omitnan');
        q = bic(:).^2; 
        Y{i,"ent2"} = -sum(q(:).*log2(q(:)),'omitnan');
        r = bic(:).^3; 
        Y{i,"ent3"} = -sum(r(:).*log2(r(:)),'omitnan');
    end

    % Copy sleep stage Annotations;
    Y.Annotations = X.Annotations;
end