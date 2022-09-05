% =======================================================
% Author Christodoulos Michaelides
% Date: August 25th, 2022
% -------------------------------------------------------
%
% Function Description: 
% This function can be used to extract features from the
% bicoherence estimations of EEG signals such as:
%   1) bicoherence entropy
%   2) bicoherence squared-entropy
%   3) average bicoherence
% -------------------------------------------------------
%
% Arguments List: (X, f)
% X: (table) a table with two columns. The first should
% contain the bicoherence estimations and the
% second should contain the sleep stage Annotations.
% You should use bicEEG() to obtain this table 
%
% f: (1D array of floats) the frequency axis of the
% bicoherence matrices. You should use bicEEG to 
% obtain this vector.
% -------------------------------------------------------
%
% Return List: (delta, theta, alpha, beta)
%
% TODO: add text
% =======================================================

function [Y] = bicoherFeatures(X, f)
    % Number of EEG segments (30sec epochs)
    N = size(X,1);

    % Initialize empty tables to store
    % features from the bicoherence matrices.
    types = ["double" "double" "double" "string"];
    names = ["ent1" "ent2" "avg" "Annotations"];
    sz = [N numel(types)];

    Y = table('Size',sz,'VariableTypes',types,'VariableNames',names);

    % binary mask to separate the primary 
    % region of the bicoherence matrix
    hex = (f >= 0) & (f' <= f) & (f + f' <= f(end));

    for i = 1:1:N
        % Extract the entire bicoherence
        % matrix for a 30sec EEG segment
        bic = cell2mat(X{i,1});
        p = bic(hex); q = p.^2;

        % Extract features from the
        % bicoherence matrix
        Y{i,"avg"}  = + mean(p);
        Y{i,"ent1"} = - mean(p .* log2(p),'omitnan');
        Y{i,"ent2"} = + log2(1 - 7.0 * mean(q .* log2(q),'omitnan'));
    end

    % Copy sleep stage annotations
    Y.Annotations = X.Annotations;
end