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
%    segment into delta, theta, alpha and beta waves.
%
% 4) Estimate the variance, skewness and kurtosis 
%    of each wave-type for every 30sec segment.
%
% 5) Plot the estimates with scatter plots. 
%    In total you should have 3 different scatter plots
%    with 3 axes each. One scatter plot for every order
%    of statistics (variance, skewness, kurtosis) and 
%    3 axes on each scatter plot (delta, theta and
%    alpha waves. Beta waves are discarded because
%    they are associated with any sleep stage)
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

channel = 1;   % EEG channel to analyse

figure(1); hold on; grid on;    % scatter plot for variance estimates
figure(2); hold on; grid on;    % scatter plot for skewness estimates
figure(3); hold on; grid on;    % scatter plot for kurtosis estimates

sz = 0.5; % size of observations on scatter plots

% Run the same analysis for 15 patients
for i = 1:1:15
    % Progress status
    fprintf('Processing file %d\n\n', i);

    % Make sure that the input files exist
    if ~isfile(sprintf("SN%03d.edf",i))
        continue;
    end

    % Load recordings from EDF files
    Z = loadEDF(i);

    % Decompose the EEG to delta, theta, alpha and beta waves
    % and estimate 2nd 3rd and 4th order statistics for every scale
    [delta, theta, alpha, beta] = mraEEG(Z,channel);

    % Update variance scatter plot
    figure(1);
    s = Z.Annotations == "Sleep stage N1";
    scatter3(delta.var(s), theta.var(s), alpha.var(s), sz, 'r');
    s = Z.Annotations == "Sleep stage N2";
    scatter3(delta.var(s), theta.var(s), alpha.var(s), sz, 'g');
    s = Z.Annotations == "Sleep stage N3";
    scatter3(delta.var(s), theta.var(s), alpha.var(s), sz, 'b');
    xlabel('delta'); ylabel('theta'); zlabel('alpha');
    legend('N1', 'N2', 'N3');
    title('scatter plot for variance estimates');

    % Update skewness scatter plot
    figure(2);
    s = Z.Annotations == "Sleep stage N1";
    scatter3(delta.skw(s), theta.skw(s), alpha.skw(s), sz, 'r');
    s = Z.Annotations == "Sleep stage N2";
    scatter3(delta.skw(s), theta.skw(s), alpha.skw(s), sz, 'g');
    s = Z.Annotations == "Sleep stage N3";
    scatter3(delta.skw(s), theta.skw(s), alpha.skw(s), sz, 'b');
    xlabel('delta'); ylabel('theta'); zlabel('alpha');
    legend('N1', 'N2', 'N3');
    title('scatter plot for skewness estimates');

    % Update kurtosis scatter plot
    figure(3);
    s = Z.Annotations == "Sleep stage N1";
    scatter3(delta.krt(s), theta.krt(s), alpha.krt(s), sz, 'r');
    s = Z.Annotations == "Sleep stage N2";
    scatter3(delta.krt(s), theta.krt(s), alpha.krt(s), sz, 'g');
    s = Z.Annotations == "Sleep stage N3";
    scatter3(delta.krt(s), theta.krt(s), alpha.krt(s), sz, 'b');
    xlabel('delta'); ylabel('theta'); zlabel('alpha');
    legend('N1', 'N2', 'N3');
    title('scatter plot for skewness estimates');
end
