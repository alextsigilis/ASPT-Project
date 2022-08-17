% ==========================================================
%
% Author: Christodoulos Michaelides
% Date: August 12th, 2022 
% ----------------------------------------------------------
% 
% Function Description: 
% This function estimates the bispectrum of an EEG signal.
% The bispectrum is estimated for every 30sec segment of 
% the EEG recording by using the direct method.
% ----------------------------------------------------------
%
% Arguments List: (X, fs, fc, K, idx)
%
% X: (table) The EEG recordings represented as a table.
% 
% Every column contains a different EEG channel. Every row of
% the table contains a 30sec segment for each EEG channel. 
% 
% You should use loadEDF() to obtain this table. You can also
% concatenate (vertically) multiple tables obtained by
% loadEDF and pass them as an argument to this function.
%
% fs: (integer) The sampling frequency of the EEG recordings
% in Hz (256Hz in our case)
%
% fc: (integer) The upper bound of the frequency axis. 
% Make sure that the Nyquist Criterion is satisfied. 
% That is: 2 * fc < fs.
%
% K: (integer) Number of partitions for every 30sec window.
% An estimation of the bispectrum is calculated for every 
% partition. The final result is obtained by averaging 
% the previous estimations. Make sure that K is not too large.
% Every partition should contain at least 5 samples, otherwise
% an exception is raised warning you against poor resolution on 
% the frequency axis.
%
% channel: (integer) An integer between 1 and 4 used as an 
% index to select an EEG channel.
% ----------------------------------------------------------
%
% Return List: (bis, freq)
%
% bis: (table) a table which contains the bispectrum estimates
% for every 30sec segment and the sleep stage labels.
%
% freq: (1D array of floats) The frequency axis of the 
% estimated bispectra. The upper bound of the axis is determined
% by the fc argument. The resolution of the axis is inversely
% proportional to the selected value of K.
% ----------------------------------------------------------
%
% A brief note on big-O time complexity and execution time:
%
% Let us assume that we are given a time-series of 
% N elements. Estimating the bispectrum with the direct
% method involves the following steps:
%
% 1) Splitting the sequence into K non-overlapping segments 
%    of length M (N = K * M)
%
% 2) Applying the FFT algorithm for each segment
%    F(i) = FFT(x) 0 <= i < M
%
% 3) Estimating the bispectrum of each segment using the 
%    following formula: B(i,j) = F(i)F(j)F*(i+j)
%    0 <= i,j < M and 0 <= i+j < M
% 
% 4) Obtaining the final estimation by averaging the
%    partial bispectra of the previous step.
%
% Step 1) requires O(N) time.
% Step 2) requires O(KMlogM) = O(NlogN - NlogK) time.
% Step 3) requires O(K*M^2) = O(N^2/K) time.
% Step 4) requires O(K*M^2) = O(N^2/K) time.
%
% As we can see, the time complexity is in the order 
% of O(N^2) assuming a constant value for K (say K=1).
%
% In our case, the sampling frequency is 
% fs = 256Hz and N = 30sec * 256Hz = 7680 samples.
% Unfortunately, this input size translates to
% tens of millions of floating point operations 
% for estimating the bispectrum of a single channel 
% for a single 30sec time segment. 
% What's even worse is that most of those operations
% are unnecessary, since they are spent estimating 
% the bispectrum for frequencies greater than 32Hz
% We are only interested in frequencies between
% 0Hz and 32Hz. 
% Therefore, it is important that you set the value 
% of fc as low as possible. Estimating the bispectrum
% with fc = 128Hz is approximately 16 times slower 
% compared to fc = 32Hz, because of the extra 
% computational overhead introduced in steps 3) 
% and 4) of the algorithm.
%
% The recommended parameters are:
% fc = 32Hz
% fs = 256Hz
% K = 10 segments
% Assuming N = 30 * 256 = 7680 samples, those 
% values give a 0.3Hz resolution in the 
% frequency axis of the bispectrum.

function [bis, freq] = bispecEEG(X,fs,fc,K,channel)
    % L: (integer) number of labeled sleep segments
    L = size(X,1);

    if L < 1 error("X is empty\n"); end
    if K < 1 error("K must be positive\n"); end
    if fs <= 0 error("fs must be positive\n"); end
    if fc <= 0 error("fc must be positive\n"); end
    if channel < 1 || channel > 4 error("Invalid EEG channel\n"); end

    % N: (integer) number of samples per labeled segment
    % M: (integer) number of samples after partitioning
    % P: (integer) number of non-negative frequencies for an M-sample FFT
    N = numel(cell2mat(X{1,channel}));
    M = floor(N/K);
    P = ceil(M/2);

    if M < 5 error("Insufficient frequency resolution\n"); end
    if 2 * fc > fs error("Nyquist criterion is violated\n"); end
    if fc * (P - 1) < fs error("fc is almost zero\n"); end

    % ind: (array of integers) array of indices for frequency axis
    % s:   (integer) upper bound of frequency axis
    ind = 1:1:(P-1);
    ind = ind(fs * ind <= 2 * fc * (P - 1)); 
    ind = ind(:);
    s = ind(end);

    if s > P - 1 error("Nyquist criterion is violated\n"); end
    if s < 2 error("fc is almost zero\n"); end

    % idx1,idx2,idx3: (arrays of integers) arrays of
    % linear indices for accumulating triple products.
    idx1 = nan;
    idx2 = nan;
    idx3 = nan;

    types = ["array" "string"];
    names = [X.Properties.VariableNames(channel) "Annotations"];

    bis = table(                ...
        'Size', [L 2],          ...
        'VariableTypes',types,  ...
        'VariableNames',names);
    
    % Copy the sleep stage annotations
    bis.Annotations = X.Annotations;

    for i = 1:1:L
        % Extract a 30sec EEG segment
        x = cell2mat(X{i,idx}); x = x(:);

        % Split x into K partitions
        x = reshape(x, [M 1 K]);

        % Subtract the mean from every partition
        x = x - mean(x,1);

        % Apply the hanning window
        x = x .* hanning(M);

        % Estimate the FFT of every partition
        x = fft(x,1);

        % Drop the DC term of every FFT
        x = x(2:end,:,:);

        % Drop frequencies higher than
        % fc and negative frequencies
        x = x(ind,:,:);

        % Accumulate triple products for bispectra
        b = x(idx1) .* x(idx2) .* x(idx3);

        % The final estimation is obtained by averaging
        % the bispectra of every partition
        bis(i,1) = sum(b,1) / K;
    end

    % Scale the indices by a factor of 
    % fs/(2*(P-1)) to obtain the frequency axis
    freq = 0.5 * fs * ind / (P - 1); 
end