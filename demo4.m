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
% of each frequency scale for every 30sec segment.
%
% 5) Plot the estimates into 3D scatter plots. 
%
% 6) Use different colors to represent each sleep stage
%
% 7) Repeat the same process for 10-20 patients in the
% original dataset
%
% 8) See if there are any high density clusters in those
% scatter plots. (Spoiler: There are)
% =======================================================

channel = 1;   % EEG channel to analyse

figure(1); hold on; grid on;    % scatter plot for delta waves
figure(2); hold on; grid on;    % scatter plot for theta waves
figure(3); hold on; grid on;    % scatter plot for alpha waves
figure(4); hold on; grid on;    % scatter plot for beta  waves

sz = 0.5; % size of observations on scatter plots

% Run the same analysis for 20 patients
for i = 1:1:20
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

    % Update scatter plot for delta wave statistics
    figure(1);
    s = delta.Annotations == "Sleep stage W";
    scatter3(delta.var(s), delta.skw(s), delta.krt(s), sz, 'r');
    s = delta.Annotations == "Sleep stage N3";
    scatter3(delta.var(s), delta.skw(s), delta.krt(s), sz, 'g');
    s = delta.Annotations == "Sleep stage R";
    scatter3(delta.var(s), delta.skw(s), delta.krt(s), sz, 'b');
    xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
    legend('W', 'N3', 'R');
    title('scatter plot for delta wave statistics');

    % Update scatter plot for theta wave statistics
    figure(2);
    s = theta.Annotations == "Sleep stage W";
    scatter3(theta.var(s), theta.skw(s), theta.krt(s), sz, 'r');
    s = theta.Annotations == "Sleep stage N3";
    scatter3(theta.var(s), theta.skw(s), theta.krt(s), sz, 'g');
    s = theta.Annotations == "Sleep stage R";
    scatter3(theta.var(s), theta.skw(s), theta.krt(s), sz, 'b');
    xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
    legend('W', 'N3', 'R');
    title('scatter plot for theta wave statistics');

    % Update scatter plot for alpha wave statistics
    figure(3);
    s = alpha.Annotations == "Sleep stage W";
    scatter3(alpha.var(s), alpha.skw(s), alpha.krt(s), sz, 'r');
    s = alpha.Annotations == "Sleep stage N3";
    scatter3(alpha.var(s), alpha.skw(s), alpha.krt(s), sz, 'g');
    s = alpha.Annotations == "Sleep stage R";
    scatter3(alpha.var(s), alpha.skw(s), alpha.krt(s), sz, 'b');
    xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
    legend('W', 'N3', 'R');
    title('scatter plot for alpha wave statistics');

    % Update scatter plot for beta wave statistics
    figure(4);
    s = beta.Annotations == "Sleep stage W";
    scatter3(beta.var(s), beta.skw(s), beta.krt(s), sz, 'r');
    s = beta.Annotations == "Sleep stage N3";
    scatter3(beta.var(s), beta.skw(s), beta.krt(s), sz, 'g');
    s = beta.Annotations == "Sleep stage R";
    scatter3(beta.var(s), beta.skw(s), beta.krt(s), sz, 'b');
    xlabel('variance'); ylabel('skewness'); zlabel('kurtosis');
    legend('W', 'N3', 'R');
    title('scatter plot for beta wave statistics');
end