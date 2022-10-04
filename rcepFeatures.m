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
%   => alpha1 (table): cepstral features for alpha waves (impulse response)
%   => alpha2 (table): cepstral features for theta waves (input sequence)
%   => beta1  (table): cepstral features for alpha waves (impulse response)
%   => beta2  (table): cepstral features for beta waves  (input sequence)
% ================================================================

function [alpha1, beta1, alpha2, beta2] = rcepFeatures(Z, channel)
    N = size(Z,1);

    % Tables to store features
    sz1    = [N 4]; 
    types1 = ["single" "single" "single" "string"];
    names1 = ["std" "skw" "krt" "Annotations"];

    sz2    = [N 4];
    types2 = ["single" "single" "single" "string"];
    names2 = ["mean" "std" "zcr" "Annotations"];

    alpha1 = table('Size',sz1,'VariableTypes',types1,'VariableNames',names1);
    alpha2 = table('Size',sz2,'VariableTypes',types2,'VariableNames',names2);

    beta1 = table('Size',sz1,'VariableTypes',types1,'VariableNames',names1);
    beta2 = table('Size',sz2,'VariableTypes',types2,'VariableNames',names2);

    for n = 1:1:N
        % Extract 30sec EEG segment
        epsilon = 1e-6;
        z = cell2mat(Z{n,channel});
        z = (z - mean(z)) ./ (std(z) + epsilon);

        % MRA decomposition
        z = modwt(z,"db6",5);
        z = modwtmra(z,"db6");

        % Extract alpha and beta waves
        win = hamming(7680); win = win(:);                
        a = z(3,:); a = a(:) .* win;                     % alpha waves
        b = z(4,:); b = b(:) .* win;                     % beta waves

        % Estimate cepstral coefficients for every frequency band
        a = real(ifft(log(epsilon + abs(fft(a)))));
        b = real(ifft(log(epsilon + abs(fft(b)))));

        % Split cepstral sequences (impulse response, input sequence)
        t = linspace(0,30,7680);
        a1 = a(0 <= t & t <= 2);
        b1 = b(0 <= t & t <= 2);
        a2 = a(2 <= t & t <= 15);
        b2 = b(2 <= t & t <= 15);

        % Feature extraction 
        alpha1{n,"std"} = std(a1);
        alpha1{n,"skw"} = skewness(a1);
        alpha1{n,"krt"} = kurtosis(a1);

        alpha2{n,"mean"} = mean(abs(a2));
        alpha2{n,"std"}  = std(a2);
        alpha2{n,"zcr"}  = zerocrossrate(detrend(a2));

        beta1{n,"std"} = std(b1);
        beta1{n,"skw"} = skewness(b1);
        beta1{n,"krt"} = kurtosis(b1);

        beta2{n,"mean"} = mean(abs(b2));
        beta2{n,"std"}  = std(b2);
        beta2{n,"zcr"}  = zerocrossrate(detrend(b2));

        % Copy sleep stage Annotations
        alpha1{n,"Annotations"} = Z{n,"Annotations"};
        alpha2{n,"Annotations"} = Z{n,"Annotations"};
        
        beta1{n,"Annotations"} = Z{n,"Annotations"};
        beta2{n,"Annotations"} = Z{n,"Annotations"};
    end
end