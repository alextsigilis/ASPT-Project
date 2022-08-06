% =======================================================
% Author: Christodoulos Michaelides
% Date: August 5th, 2022
% -------------------------------------------------------
%
% Objectives:
% 1) Choose one out of the 4 EEG channels.
%
% 2) Split the recording of the EEG into 30sec segments.
%
% 3) Use Multi Resolution Analysis to decompose every 
% segment into delta, theta, alpha and beta waves.
%
% 4) Estimate the variance, skewness and kurtosis 
%    of each wave-type for every 30sec segment.
%
% 5) Plot the estimates with scatter plots. 
%
% 6) Use different colors to represent each sleep stage
%
% 7) Repeat the same process for 10-20 patients in the
% original dataset
%
% 8) See if there are any high-density clusters in those
% scatter plots. (Spoiler: There are)
% =======================================================

clear all; close all; clc; % Reset your workspace

channel = 1;               % EEG channel to analyse
num_of_patients = 25;      % number of patients 

% -------------------------------------------------------
% Do not change anything below that point
% -------------------------------------------------------

sz = 0.25; % size of observations on scatter plots

% List of valid sleep stages
stages = [
    "Sleep stage W",  ... 
    "Sleep stage N1", ...
    "Sleep stage N2", ...
    "Sleep stage N3", ...
    "Sleep stage R"];

num_of_stages = length(stages);

% an incremental index for every plot
idx = 1;

% Initialize all scatter plots
% (axes titles, legends etc ...)
for i = 1:num_of_stages
    for j = (i+1):num_of_stages
        % delta waves
        figure(idx); hold on; grid on;
        title(sprintf("Delta waves: %s vs %s", stages(i), stages(j)));
        xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
        idx = idx + 1;

        % theta waves 
        figure(idx); hold on; grid on;
        title(sprintf("Theta waves: %s vs %s", stages(i), stages(j)));
        xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
        idx = idx + 1;

        % alpha waves
        figure(idx); hold on; grid on;
        title(sprintf('Alpha waves: %s vs %s', stages(i), stages(j)));
        xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
        idx = idx + 1;

        % beta waves
        figure(idx); hold on; grid on;
        title(sprintf('Beta waves: %s vs %s', stages(i), stages(j)));
        xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
        idx = idx + 1;
    end
end

for n = 1:1:num_of_patients
    % Progress status
    fprintf('Processing file %d out of %d\n\n', n, num_of_patients);

    % Make sure that the input files exist
    if ~isfile(sprintf("SN%03d.edf",n))
        continue;
    end

    % Load recordings from EDF files
    Z = loadEDF(n);

    % Decompose the EEG to delta, theta, alpha and beta waves
    % and estimate 2nd, 3rd and 4th order statistics for every scale
    [delta, theta, alpha, beta] = mraEEG(Z,channel);

    % An incremental index for every scatter plot
    idx = 1;
    
    % Update all scatter plots with new data
    for i = 1:num_of_stages
        for j = (i+1):num_of_stages
            % scatter plot for delta waves
            figure(idx);
            s = Z.Annotations == stages(i);
            scatter3(delta.var(s), delta.skw(s), delta.krt(s), sz, 'r');
            s = Z.Annotations == stages(j);
            scatter3(delta.var(s), delta.skw(s), delta.krt(s), sz, 'g');
            idx = idx + 1;

            % scatter plot for theta waves
            figure(idx);
            s = Z.Annotations == stages(i);
            scatter3(theta.var(s), theta.skw(s), theta.krt(s), sz, 'r');
            s = Z.Annotations == stages(j);
            scatter3(theta.var(s), theta.skw(s), theta.krt(s), sz, 'g');
            idx = idx + 1;

            % scatter plot for alpha waves
            figure(idx);
            s = Z.Annotations == stages(i);
            scatter3(alpha.var(s), alpha.skw(s), alpha.krt(s), sz, 'r');
            s = Z.Annotations == stages(j);
            scatter3(alpha.var(s), alpha.skw(s), alpha.krt(s), sz, 'g');
            idx = idx + 1;

            % scatter plot for beta waves
            figure(idx);
            s = Z.Annotations == stages(i);
            scatter3(beta.var(s), beta.skw(s), beta.krt(s), sz, 'r');
            s = Z.Annotations == stages(j);
            scatter3(beta.var(s), beta.skw(s), beta.krt(s), sz, 'g');
            idx = idx + 1;
        end
    end
end

% incremental index for scatter plots
idx = 1;

% add a legend on every scatter plot
for i = 1:num_of_stages
    for j = (i+1):num_of_stages
        % legend for delta waves
        figure(idx);
        legend(stages(i), stages(j));
        idx = idx + 1;

        % legend for theta waves
        figure(idx);
        legend(stages(i), stages(j));
        idx = idx + 1;
        
        % legend for alpha waves
        figure(idx);
        legend(stages(i), stages(j));
        idx = idx + 1;
        
        % legend for beta waves
        figure(idx);
        legend(stages(i), stages(j));
        idx = idx + 1;
    end
end