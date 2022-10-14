% Orfanos Dimitrios, 9579

% Scatter plots for the features for the ECG

clear all; close all; clc;
tic
channel = 8;   % ECG channel to analyse
fs = 256;
gr = 0;

figure(1); hold on; grid on;     % scatter plot for the W - N1 sleep stages
figure(2); hold on; grid on;     % scatter plot for the W - N2 sleep stages
figure(3); hold on; grid on;     % scatter plot for the W - N3 sleep stages
figure(4); hold on; grid on;     % scatter plot for the W - R sleep stages
figure(5); hold on; grid on;     % scatter plot for the N1 - N2 sleep stages
figure(6); hold on; grid on;     % scatter plot for the N1 - N3 sleep stages
figure(7); hold on; grid on;     % scatter plot for the N1 - R sleep stages
figure(8); hold on; grid on;     % scatter plot for the N2 - N3 sleep stages
figure(9); hold on; grid on;     % scatter plot for the N2 - R sleep stages
figure(10); hold on; grid on;    % scatter plot for the N3 - R sleep stages
figure(11); hold on; grid on;    %%%%%%%



sz = 0.3; % size of observations on scatter plots

% Run the same analysis for 10 patients
for i = 1:11
    % Progress status
    fprintf('Processing file %d\n\n', i);

    % Make sure that the input files exist
    if ~isfile(sprintf("SN%03d.edf",i))
        continue;
    end

    % Load recordings from EDF files
    Z = loadEDF(i);
    
    % ecg_feat = hrhrv_pan_tompkin(i, fs, gr)
    ecg_feat = featuresECG(i, fs, gr);

    %% Update scatter plot for the W - N1 sleep stages
    figure(1);
    s = ecg_feat.Annotations == "Sleep stage W";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage N1";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('W', 'N1');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the W - N2 sleep stages
    figure(2);
    s = ecg_feat.Annotations == "Sleep stage W";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage N2";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('W', 'N2');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the W - N3 sleep stages
    figure(3);
    s = ecg_feat.Annotations == "Sleep stage W";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage N3";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('W', 'N3');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the W - R sleep stages
    figure(4);
    s = ecg_feat.Annotations == "Sleep stage W";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage R";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('W', 'R');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the N1 - N2 sleep stages
    figure(5);
    s = ecg_feat.Annotations == "Sleep stage N1";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage N2";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('N1', 'N2');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the N1 - N3 sleep stages
    figure(6);
    s = ecg_feat.Annotations == "Sleep stage N1";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage N3";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('N1', 'N3');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the N1 - R sleep stages
    figure(7);
    s = ecg_feat.Annotations == "Sleep stage N1";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage R";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('N1', 'R');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the N2 - N3 sleep stages
    figure(8);
    s = ecg_feat.Annotations == "Sleep stage N2";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage N3";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('N2', 'N3');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the N2 - R sleep stages
    figure(9);
    s = ecg_feat.Annotations == "Sleep stage N2";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage R";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('N2', 'R');
    title('Scatter plot for ECG Features');
    
    %% Update scatter plot for the N3 - R sleep stages
    figure(10);
    s = ecg_feat.Annotations == "Sleep stage N3";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'r');
    s = ecg_feat.Annotations == "Sleep stage R";
    scatter(ecg_feat.heartrate(s), ecg_feat.HeartRateVariability(s), sz, 'b');
    xlabel('Heart Rate'); ylabel('Heart Rate Variability');
    legend('N3', 'R');
    title('Scatter plot for ECG Features');
    
    figure(11);
end
toc