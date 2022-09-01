function [bic, freq] = bicEEGfast(X, K, fs, fc, channel, method)

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
    % len: length of truncated FFT
    % win: hanning window for FFT
    % tri: array of indices for accumulating triple products
    % seg: array of indices for partitioning the EEG records
    idx = 1:M; win = hanning(M);
    idx = idx(-fc * M <= freq * fs & freq * fs <= fc * M);
    len = numel(idx);
    tri = hankel([1:len],[len,1:len-1]);                                    % TODO: repmat, sub2ind, ind2sub?
    seg = [1:M]' + M*[0:K-1];

    % Hexagonal mask to remove artifacts
    % outside the hexagonal symmetry regions
    lo = -(len-1)/2; hi = +(len-1)/2; u = lo:1:hi;
    hex = ones(len,1) * u;
    hex = abs(hex) + abs(hex') + abs(hex+hex');
    hex = (hex < len);
    hex = repmat(hex, [1 1 K]);

    % epsilon: (float) a small positive constant to
    % ensure numerical stability when performing
    % floating point divisions
    epsilon = 1e-5;

    % Estimate the bispectrum of every 30sec EEG record
    for i = 1:1:N
        % Extract a 30sec EEG record
        x = cell2mat(X{i,channel});

        % Split, standardize and apply a window function
        y = x(seg);
        y = (y - mean(y,1)) ./ (std(y,0,1) + epsilon);           
        y = y .* win;

        % estimate the FFT of every segment
        % and discard unnecessary frequencies
        Y = fft(y,1) / M;
        Y = fftshift(Y,1);
        Y = Y(idx,:);
        Y = ifftshift(Y,1);
        CY = conj(Y);

        % Reshape the FFT matrix appropriately in order
        % to vectorize the following calculations
        Y1 = reshape(Y, [M 1 K]);                                           
        Y2 = reshape(Y, [1 M K]);

        % Estimate the bispectrum (b) and the 
        % normalization coefficients (Y12, Y3)
        Y12 = Y1 .* Y2; 
        Y3  = CY(tri); 
        b   = Y12 .* CY;

        % Normalize the bispectrum matrix to obtain
        % the bicoherence.
        Y12 = abs(Y12) .^ 2; 
        Y12 = sum(Y12,3);

        Y3 = abs(Y3) .^ 2;
        Y3 = sum(Y3,3);

        b = sum(b,3);
        b = abs(b) .^ 2;

        b = b ./ (Y12 .* Y3 + epsilon);

        % Shift the elements of the bispectrum matrix,
        % discard any elements outside the symmetry 
        % regions and save the final result.
        b = fftshift(b);
        bic{i,1} = {b .* hex};
    end

    % ---------------------------------------------------------------
    % Truncate the frequency axis and rescale according to fs 
    % ---------------------------------------------------------------
    freq = freq(idx);
    freq = freq * fs / M;
end