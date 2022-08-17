% ====================================================================
% Author: Christodoulos Michaelides
% Date: August 5th, 2022
% --------------------------------------------------------------------
%
% Objectives:
%
% 1) Choose a number of patients from the original dataset and an
% EEG channel
%
% 2) Decompose the EEG recordings of those patients into different 
% frequency scales by applying the Discrete Wavelet Transform
%
% 3) Estimate the variance, skewness and (excess) kurtosis of the
% DWT coefficients for every frequency scale
%
% 4) Use 3D scatterplots to visualize the distributions of those
% quantities
%
% ====================================================================

% ====================================================================
% 0) Reset your Workspace
% ====================================================================
clear all; close all; clc;

% ====================================================================
% 1) Script parameters (Choose any value you want to experiment with)
% ====================================================================

% n1: (integer 1-154) ID of first patient 
% n2: (integer 1-154) ID of last patient
% channel: (integer 1-4) ID of EEG channel 
% sz: (float) marker size for scatter plots
% path: save path for scatter plots
n1 = 1; n2 = 20;
channel = 4;
sz = 0.35;
path = "C:\Users\USER\Desktop\figures4\";

% --------------------------------------------------------------------
% Do not change anything below that point, unless you know what you 
% you are doing.
% --------------------------------------------------------------------

% names: (array) array of column names for table X
% types: (array) array of datatypes for every colunn of X
names = {'var'     'skw'     'krt'     'Annotations'};
types = {'double'  'double'  'double'  'string'};

delta = table('Size', [0 4], 'VariableTypes', types, 'VariableNames', names);
theta = table('Size', [0 4], 'VariableTypes', types, 'VariableNames', names);
alpha = table('Size', [0 4], 'VariableTypes', types, 'VariableNames', names);
beta  = table('Size', [0 4], 'VariableTypes', types, 'VariableNames', names);

% ====================================================================
% 2) Decompose the input signals by using DWT. Extract statistical 
%    features such as the standard deviation, skewness and kurtosis
%    of the DWT coefficients for every frequency scale.
% ====================================================================

for i = n1:1:n2
    % Progress status
    fprintf("Processing patient %d ...\n",i);

    % Make sure the input file exists
    % before attempting to open it
    if ~isfile(sprintf("SN%03d.edf",i))
        fprintf("Patient %d does not exist\n\n",i);
        continue;
    end

    % a) Load the EEG recordings and sleep stage labels
    % b) Estimate DWT coefficients for every frequency scale
    % c) Extract features such as std, skewness and kurtosis
    %    from the DWT coefficients
    % d) save the results and move on to the next patient
    
    fprintf("Loading EEG recordings from disk ...");
    Z = loadEDF(i);
    fprintf("Done\n");

    fprintf("Extracting DWT coefficients ...");
    Z = mraEEG(Z,channel);
    fprintf("Done\n");

    fprintf("Extracting features from DWT coefficients ...");
    [d, u, a, b] = statEEG(Z);
    fprintf("Done\n\n");
    
    delta = vertcat(delta, d);
    theta = vertcat(theta, u);
    alpha = vertcat(alpha, a);
    beta  = vertcat(beta,  b); 
end

% ====================================================================
% 3) Use 3-dimensional scatter plots to visualize the distributions 
%    of the extracted features. See if any of the following happen 
%    by visually inspecting the scatter plots:
%    a) Formation of high density clusters in 3D space
%    b) Separation of said clusters based on the sleep stage
% ====================================================================

% stages: (1D-array) array of valid sleep stage Annotations
stages = [                  ...
    "Sleep stage W",        ...
    "Sleep stage R",        ...
    "Sleep stage N1",       ...
    "Sleep stage N2",       ...
    "Sleep stage N3"];

% Incremental index for scatter plots
idx = 1;

for i = 1:1:numel(stages)
    for j = (i+1):1:numel(stages)
        % scatter plot for delta waves
        f = figure(idx); hold on; grid on;
        mask1 = delta.Annotations == stages(i);
        mask2 = delta.Annotations == stages(j);
        var1 = delta.var(mask1);
        var2 = delta.var(mask2);
        skw1 = delta.skw(mask1);
        skw2 = delta.skw(mask2);
        krt1 = delta.krt(mask1);
        krt2 = delta.krt(mask2);
        scatter3(var1, skw1, krt1, sz, [0.0 0.0 1.0]);
        scatter3(var2, skw2, krt2, sz, [1.0 0.6 0.0]);
        xlabel('Variance of DWT coefficients');
        ylabel('Skewness of DWT coefficients');
        zlabel('Kurtosis of DWT coefficients');
        title(sprintf("Delta waves %s vs %s", stages(i), stages(j)));
        legend(stages(i), stages(j));
        saveas(f,sprintf("%s%02d.fig",path,idx));
        idx = idx + 1;

        % scatter plot for theta waves
        f = figure(idx); hold on; grid on;
        mask1 = theta.Annotations == stages(i);
        mask2 = theta.Annotations == stages(j);
        var1 = theta.var(mask1);
        var2 = theta.var(mask2);
        skw1 = theta.skw(mask1);
        skw2 = theta.skw(mask2);
        krt1 = theta.krt(mask1);
        krt2 = theta.krt(mask2);
        scatter3(var1, skw1, krt1, sz, [0.0 0.0 1.0]);
        scatter3(var2, skw2, krt2, sz, [1.0 0.6 0.0]);
        xlabel('Variance of DWT coefficients');
        ylabel('Skewness of DWT coefficients');
        zlabel('Kurtosis of DWT coefficients');
        title(sprintf("Theta waves %s vs %s", stages(i), stages(j)));
        legend(stages(i), stages(j));
        saveas(f,sprintf("%s%02d.fig",path,idx));
        idx = idx + 1;

        % scatter plot for alpha waves
        f = figure(idx); hold on; grid on;
        mask1 = alpha.Annotations == stages(i);
        mask2 = alpha.Annotations == stages(j);
        var1 = alpha.var(mask1);
        var2 = alpha.var(mask2);
        skw1 = alpha.skw(mask1);
        skw2 = alpha.skw(mask2);
        krt1 = alpha.krt(mask1);
        krt2 = alpha.krt(mask2);
        scatter3(var1, skw1, krt1, sz, [0.0 0.0 1.0]);
        scatter3(var2, skw2, krt2, sz, [1.0 0.6 0.0]);
        xlabel('Variance of DWT coefficients');
        ylabel('Skewness of DWT coefficients');
        zlabel('Kurtosis of DWT coefficients');
        title(sprintf("Alpha waves %s vs %s", stages(i), stages(j)));
        legend(stages(i), stages(j));
        saveas(f,sprintf("%s%02d.fig",path,idx));
        idx = idx + 1;

        % scatter plot for beta  waves
        f = figure(idx); hold on; grid on;
        mask1 = beta.Annotations == stages(i);
        mask2 = beta.Annotations == stages(j);
        var1 = beta.var(mask1);
        var2 = beta.var(mask2);
        skw1 = beta.skw(mask1);
        skw2 = beta.skw(mask2);
        krt1 = beta.krt(mask1);
        krt2 = beta.krt(mask2);
        scatter3(var1, skw1, krt1, sz, [0.0 0.0 1.0]);
        scatter3(var2, skw2, krt2, sz, [1.0 0.6 0.0]);
        xlabel('Variance of DWT coefficients');
        ylabel('Skewness of DWT coefficients');
        zlabel('Kurtosis of DWT coefficients');
        title(sprintf("Beta waves %s vs %s", stages(i), stages(j)));
        legend(stages(i), stages(j));
        saveas(f,sprintf("%s%02d.fig",path,idx));
        idx = idx + 1;
    end
end