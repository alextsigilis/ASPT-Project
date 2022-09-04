% ===================================================================
% Author: Christodoulos Michaelides
% Date: August 17th, 2022
% -------------------------------------------------------------------
%
% Function Description:
% This function estimates the bicoherence index of EEG recordings. 
% By locating the peaks of the bicoherence, we can detect quadratic
% phase coupling between different frequencies. The bicoherence is 
% essentially a normalized bispectrum. Many methods have been used
% to normalize the bispectrum. This function uses the one proposed 
% by Nagashima, 2006:
% 
%                 E{|F(f1)F(f2)F'(f1+f2)|^2}      
% b(f1,f2) = ------------------------------------  (1)
%             E{|F(f1)F(f2)|^2} E{|F'(f1+f2)|^2}
%
% where:
%   1) E{*} denotes mathematical expectation
%   2) the apostrophe ' denotes complex conjugation
%   3) F is the FFT of an EEG segment
%   4) b is the bicoherence index
%
% The estimation of the bicoherence is obtained in the following
% manner:
%
% 1) The EEG signal is split into K partitions.
% 2) The samples are standardized by subtracting the mean and 
%    dividing with the standard deviation of every partition.
% 3) A window function (hanning window) is applied in every 
%    partition.
% 4) The bispectrum of every partition is estimated, as well as 
%    some normalization coefficients related to the power spectrum.
%    See formula (1).
% 5) The final estimation of the bicoherence is obtained by averaging
%    the bispectrum and normalization coefficients according to 
%    formula (1).
% -------------------------------------------------------------------
%
% Arguments List: (X, K, fs, fc, channel, method)
%
% X: (table) the EEG recordings and sleep stage
% Annotations. You should use loadEDF to obtain this 
% table.
%
% K: (integer) number of partitions for estimating 
% the bicoherence in a 30sec epoch. The bispectrum is 
% estimated for every partition. The bicoherence is 
% obtained by normalizing the bispectra across the entire
% 30sec epoch.
%
% fs: (float) sampling frequency of EEG recordings
% (Usually 256Hz).
%
% fc: (float) The upper limit of the frequency axis.
% Since it is rarely necessary to estimate the
% bispectrum for every pair of frequencies
% -fs/2 < f1,f2 < +fs/2, we can greatly reduce the 
% execution time and memory requirements by truncating
% the frequency axis. Always make sure that fs and 
% fc satisfy the Nyquist criterion. That is: 2*fc < fs
%
% channel: (integer) An integer between 1 and 4 used to 
% select one out of the 4 available EEG channels.
%
% method: (string) You can choose between "fancy" and "fast".
% The first method ("fancy") estimates the bicoherence matrix
% both in its primary region and the symmetric ones. This is 
% the prefered method when we are interested in obtaining 
% visually pleasing results (contour plots). The second method
% ("fast") skips the symmetric regions entirely and focuses only on
% the primary one. This is the prefered method when we are 
% interested in estimating the bicoherence for a large number of 
% patients and EEG channels. The "fast" method should be 
% approximately 4 times faster and more memory-efficient compared to 
% the "fancy" one.
% -------------------------------------------------------------------
%
% Return Variables: (bic, freq)
%
% bic: (table) A table with two columns. The first column
% stores the bicoherence estimations for every 30sec epoch.
% The second column stores the sleep stage Annotations 
%
% freq: (array of floats) A 1D array which can be used as a 
% frequency axis for the bicoherence matrix. The size and the 
% bounds of freq depend on fs, fc, K and method.
% ===================================================================

function [bic, freq] = bicEEG(X, K, fs, fc, channel, method)

    % ---------------------------------------------------------------
    % Parameter Checks
    % ---------------------------------------------------------------

    if K  <= 0 error("K must be positive "); end
    if fs <= 0 error("fs must be positive"); end
    if fc <= 0 error("fc must be positive"); end
    if channel < 1 || channel > 4 error("invalid EEG channel"); end

    if method ~= "fancy" && method ~= "fast"
        error("Invalid estimation method");
    end

    N = size(X,1);
    if N <= 0 error("X is empty"); end

    L = numel(cell2mat(X{1,channel}));
    if L <= 0 error("EEG records are empty"); end

    M = floor(L/K);
    if M < 5 error("low frequency resolution"); end
    if 2*M*fc >= (M-2)*fs error("Nyquist criterion is violated"); end
    if fc*M <= fs error("fc is almost zero"); end

    % ---------------------------------------------------------------
    % Create a table to store the bicoherence and sleep stage labels
    % ---------------------------------------------------------------

    names = [X.Properties.VariableNames{channel}, "Annotations"];
    types = ["cell", "string"];
    sz = [N numel(types)];
    bic = table('Size',sz,'VariableTypes',types,'VariableNames',names);
    bic.Annotations = X.Annotations;

    % ---------------------------------------------------------------
    % Frequency axis
    % ---------------------------------------------------------------

    if rem(M,2) == 0
        freq = [-M/2:(M/2-1)];
    elseif rem(M,2) == 1
        freq = [-(M-1)/2:(M-1)/2];
    end

    % ---------------------------------------------------------------
    % Estimation of Bispectra
    % ---------------------------------------------------------------

    % idx: array of indices for discarding unnecessary FFT components
    %
    % len: length of truncated FFT
    %
    % win: hanning window for FFT
    %
    % tri: array of indices for accumulating triple products
    %
    % seg: array of indices for partitioning the EEG records
    %
    % hex: hexagonal boolean mask to remove artifacts outside the
    % symmetry regions
    %
    % epsilon: a small positive constant to ensure numerical
    % stability when performing floating point divisions

    idx = 1:M; win = hanning(M);

    if method == "fancy"
        idx = idx(-fc * M <= freq * fs & freq * fs <= fc * M);
    elseif method == "fast"
        idx = idx(freq * fs <= fc * M & freq >= 0);
    end

    len = numel(idx);

    tri = hankel([1:len],[len,1:len-1]);
    tri = reshape(tri, [len len 1]) + reshape(len * [0:K-1], [1 1 K]);

    seg = [1:M]' + M*[0:K-1];

    if method == "fancy"
        u = (1-len)/2:1:(len-1)/2;
        u = ones(len,1) * u;
        u = abs(u) + abs(u') + abs(u + u');
        hex = u < len;
    elseif method == "fast"
        u   = 0:1:(len-1);
        u   = ones(len,1) * u;
        hex = (u' <= u) & (u >= 0) & (u + u' < len);
    end

    epsilon = 1e-5;

    % Estimate the bispectrum of every 30sec EEG record
    if method == "fancy"
        for i = 1:1:N
            % Extract a 30sec EEG record
            x = cell2mat(X{i,channel});

            % Split, standardize and apply a window function
            y = x(seg);
            y = (y - mean(y,1)) ./ (std(y,0,1) + epsilon);
            y = y .* win;

            % estimate the FFT of every segment
            % and discard unnecessary frequencies
            Y = fft(y,[],1) / M;
            Y = fftshift(Y,1);
            Y = Y(idx,:);
            Y = ifftshift(Y,1);
            CY = conj(Y); CY = CY(tri);

            % Reshape the FFT matrix appropriately in order
            % to vectorize the following calculations
            Y1 = reshape(Y, [len 1 K]);
            Y2 = reshape(Y, [1 len K]);

            % Estimate the bispectrum (b) and the
            % normalization coefficients (Y12, CY)
            Y12 = Y1 .* Y2;
            b   = Y12 .* CY;

            % Normalize the bispectrum matrix to obtain
            % the bicoherence.
            Y12 = abs(Y12) .^ 2;
            Y12 = sum(Y12, 3);

            CY = abs(CY) .^ 2;
            CY = sum(CY, 3);

            b = sum(b, 3);
            b = abs(b) .^ 2;

            b = b ./ (Y12 .* CY + epsilon);

            % Shift the elements of the bispectrum matrix,
            % discard any elements outside the symmetry
            % regions and save the final result.
            b = fftshift(b);
            bic{i,1} = {b .* hex};
        end

    elseif method == "fast"
        for i = 1:1:N
            % Extract a 30sec EEG record
            x = cell2mat(X{i,channel});
    
            % Split, standardize and apply a window function
            y = x(seg);
            y = (y - mean(y,1)) ./ (std(y,0,1) + epsilon);
            y = y .* win;
    
            % estimate the FFT of every segment
            % and discard unnecessary frequencies
            Y = fft(y,[],1) / M;
            Y = fftshift(Y,1);
            Y = Y(idx,:);
            CY = conj(Y); CY = CY(tri);
    
            % Reshape the FFT matrix appropriately in order
            % to vectorize the following calculations
            Y1 = reshape(Y, [len 1 K]);
            Y2 = reshape(Y, [1 len K]);
    
            % Estimate the bispectrum (b) and the
            % normalization coefficients (Y12, CY)
            Y12 = Y1 .* Y2;
            b   = Y12 .* CY;
    
            % Normalize the bispectrum matrix to obtain
            % the bicoherence.
            Y12 = abs(Y12) .^ 2;
            Y12 = sum(Y12, 3);
    
            CY = abs(CY) .^ 2;
            CY = sum(CY, 3);
    
            b = sum(b, 3);
            b = abs(b) .^ 2;
    
            b = b ./ (Y12 .* CY + epsilon);
    
            % Discard any elements outside the symmetry
            % regions and save the final result.
            bic{i,1} = {b .* hex};
        end
    end

    % ---------------------------------------------------------------
    % Truncate the frequency axis and rescale according to fs
    % ---------------------------------------------------------------
    freq = freq(idx);
    freq = freq * fs / M;
end