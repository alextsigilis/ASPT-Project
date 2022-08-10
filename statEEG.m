% =================================================
% Author: Christodoulos Michaelides
% Date: August 10th, 2022
% -------------------------------------------------
%
% Function Description:
% This function can be used to extract statistical
% features from the DWT domain such as the 
% standard deviation, skewness and kurtosis of the
% DWT coefficients. This process is carried out
% for every DWT scale independently. The resulting
% features can then be used by a classifier to 
% distinguish between different sleep stages
% (REM sleep, Non-REM sleep etc ...)
% -------------------------------------------------
%
% Arguments List: (X)
%
% X: (table) a table with 4 + 1 columns obtained 
% either by calling the mraEEG() function for a
% single patient or by vertically-concatenating 
% multiple such tables for different patients.
%
% The 1st column contains the DWT coefficients
% of the delta waves 
% The 2nd column contains the DWT coefficients
% of the theta waves
% The 3rd column contains the DWT coefficients 
% of the alpha waves
% The 4th column contains the DWT coefficients 
% of the beta waves
% The 5th column contains the sleep-stage
% Annotations for every 30sec segment.
%
% You should read the documentation of mraEEG.m
% for more information
% -------------------------------------------------
% 
% Return Variables: (delta, theta, alpha, beta)
%
% delta: (table) a table with 3 + 1 columns.
% The 1st column contains the standard deviation
% estimated using the DWT coefficients of the 
% delta waves.
% The 2nd column contains the skewness estimated
% using the DWT coefficients of the delta waves
% The 3rd column contains the kurtosis estimated
% using the DWT coefficients of the delta waves
% The 4th column contains the sleep-stage
% Annotations for every 30sec segment.
%
% theta: (table) same as delta but for theta waves
% instead.
%
% alpha: (table) same as delta but for alpha waves
% instead.
%
% beta: (table) same as delta but for beta waves
% instead.
%
% -------------------------------------------------
%
% This is what the output tables should look like:
%
%      var        skw        krt        Annotations    
%    ______    _________    ______    _______________
%
%    78.617     0.32108     3.8940    "Sleep stage W"
%    71.075    -0.42419     3.7714    "Sleep stage N1"
%      :           :          :              :       
%    403.52     0.44983     11.019    "Sleep stage W"
% =================================================

function [delta, theta, alpha, beta] = statEEG(X)
    % names: Array of names for every column
    names = ["var" "skw" "krt" "Annotations"];

    % varFunc: function handle to estimate standard deviation
    % skwFunc: function handle to estimate skewness
    % krtFunc: function handle to estimate (excess) kurtosis
    varFunc = @(x) std(cell2mat(x)); 
    skwFunc = @(x) skewness(cell2mat(x));
    krtFunc = @(x) kurtosis(cell2mat(x)) - 3;

    % Estimate statistics for delta waves
    delta(:,1) = rowfun(varFunc,X(:,"delta"));
    delta(:,2) = rowfun(skwFunc,X(:,"delta"));
    delta(:,3) = rowfun(krtFunc,X(:,"delta"));
    delta = addvars(delta, X.Annotations);

    % Estimate statistics for theta waves
    theta(:,1) = rowfun(varFunc,X(:,"theta"));
    theta(:,2) = rowfun(skwFunc,X(:,"theta"));
    theta(:,3) = rowfun(krtFunc,X(:,"theta"));
    theta = addvars(theta, X.Annotations);

    % Estimate statistics for alpha waves
    alpha(:,1) = rowfun(varFunc, X(:,"alpha"));
    alpha(:,2) = rowfun(skwFunc, X(:,"alpha"));
    alpha(:,3) = rowfun(krtFunc, X(:,"alpha"));
    alpha = addvars(alpha, X.Annotations);

    % Estimate statistics for beta waves
    beta(:,1) = rowfun(varFunc, X(:,"beta"));
    beta(:,2) = rowfun(skwFunc, X(:,"beta"));
    beta(:,3) = rowfun(krtFunc, X(:,"beta"));
    beta = addvars(beta, X.Annotations);

    % Add variable names to every column
    oldnames = delta.Properties.VariableNames;
    delta = renamevars(delta, oldnames, names);
    oldnames = theta.Properties.VariableNames;
    theta = renamevars(theta, oldnames, names);
    oldnames = alpha.Properties.VariableNames;
    alpha = renamevars(alpha, oldnames, names);
    oldnames = beta.Properties.VariableNames;
    beta  = renamevars(beta,  oldnames, names);
end