% ================================================================
% Authors:  Chrysa Doulou, 
%           Christodoulos Michaelides,
%           Demetrios Orphanos
%           
% Date: September 30th, 2022
% ----------------------------------------------------------------
%
% Function Description:
% This function splits EEG recordings into delta, theta, alpha 
% and beta waves. It then estimates the cepstal coefficients 
% of every frequency scale and extracts a set of features from
% every cepstral sequence. These features include:
%   => The average of cepstral coefficients
%   => The standard deviation of cepstral coefficients
%   => The skewness of cepstral coefficients
%   => The kurtosis of cepstral coefficients
%   => The zero-crossing rate of cepstral coefficients
% ----------------------------------------------------------------
%
% Arguments List:
%   => Z (table) a table which contains 30sec epochs of EEG
%   recordings. You can use loadEDF() to obtain the raw data
%   of the EEG recordings or prefilter() to apply an elliptic
%   filter to the EEG recordings first.
%
%   => channel (string): The EEG channel we are interested in.
%   For example "EEGC3_M2", "EEGF4_M1" etc...
% ----------------------------------------------------------------
%
% Return variables:
%   => delta (table): cepstral features for delta waves
%   => theta (table): cepstral features for theta waves
%   => alpha (table): cepstral features for alpha waves
%   => beta  (table): cepstral features for beta waves
% ================================================================

function [delta, theta, alpha, beta] = rcepFeatures(Z, channel)
    N = size(Z,1);

    sz = [N 6]; 
    types = ["single" "single" "single" "single" "single" "string"];
    names = ["mean" "std" "skw" "krt" "zcr" "Annotations"];

    delta = table('Size',sz,'VariableTypes',types,'VariableNames',names);
    theta = table('Size',sz,'VariableTypes',types,'VariableNames',names);
    alpha = table('Size',sz,'VariableTypes',types,'VariableNames',names);
    beta  = table('Size',sz,'VariableTypes',types,'VariableNames',names);

    for n = 1:1:N
        % Extract 30sec EEG segment
        epsilon = 1e-6;
        z = cell2mat(Z{n,channel});
        z = (z - mean(z)) ./ (std(z) + epsilon);

        % MRA decomposition
        z = modwt(z,"db6",5);
        z = modwtmra(z,"db6");

        % Extract delta, theta, alpha and beta waves
        win = hamming(7680); win = win(:);
        z1 = z(1,:); z1 = z1(:) .* win;                     % delta waves                   
        z2 = z(2,:); z2 = z2(:) .* win;                     % theta waves
        z3 = z(3,:); z3 = z3(:) .* win;                     % alpha waves
        z4 = z(4,:); z4 = z4(:) .* win;                     % beta waves

        % Estimate cepstral coefficients for every frequency band
        z1 = real(ifft(log(epsilon + abs(fft(z1)))));
        z2 = real(ifft(log(epsilon + abs(fft(z2)))));
        z3 = real(ifft(log(epsilon + abs(fft(z3)))));
        z4 = real(ifft(log(epsilon + abs(fft(z4)))));

        % Crop symmetric regions from cepstral sequences
        t = linspace(0,30,7680);
        z1 = z1(0 <= t & t <= 2);
        z2 = z2(0 <= t & t <= 2);
        z3 = z3(0 <= t & t <= 2);
        z4 = z4(0 <= t & t <= 2);
        % z1 = z1(2 <= t & t <= 15);
        % z2 = z2(2 <= t & t <= 15);
        % z3 = z3(2 <= t & t <= 15);
        % z4 = z4(2 <= t & t <= 15);

        % absolute mean of cepstral coefficients
        delta{n,"mean"} = mean(abs(z1));
        theta{n,"mean"} = mean(abs(z2));
        alpha{n,"mean"} = mean(abs(z3));
        beta{n,"mean"}  = mean(abs(z4));

        % standard deviation of cepstral coefficients
        delta{n,"std"} = std(z1);
        theta{n,"std"} = std(z2);
        alpha{n,"std"} = std(z3);
        beta{n,"std"}  = std(z4);

        % skewness of cepstral coefficients
        delta{n,"skw"} = skewness(z1);
        theta{n,"skw"} = skewness(z2);
        alpha{n,"skw"} = skewness(z3);
        beta{n,"skw"}  = skewness(z4);

        % kurtosis of cepstral coefficients (log scale)
        delta{n,"krt"} = log(1 + kurtosis(z1));
        theta{n,"krt"} = log(1 + kurtosis(z2));
        alpha{n,"krt"} = log(1 + kurtosis(z3));
        beta{n,"krt"}  = log(1 + kurtosis(z4));

        % zero crossing rate of cepstral coefficients
        delta{n,"zcr"} = zerocrossrate(detrend(z1,1));
        theta{n,"zcr"} = zerocrossrate(detrend(z2,1));
        alpha{n,"zcr"} = zerocrossrate(detrend(z3,1));
        beta{n,"zcr"}  = zerocrossrate(detrend(z4,1));

        % Copy sleep stage Annotations
        delta{n,"Annotations"} = Z{n,"Annotations"};
        theta{n,"Annotations"} = Z{n,"Annotations"};
        alpha{n,"Annotations"} = Z{n,"Annotations"};
        beta{n,"Annotations"}  = Z{n,"Annotations"};
    end
end