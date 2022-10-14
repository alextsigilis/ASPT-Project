% Orfanos Dimitrios, 9579

% Heart Rate & Heart Rate Variability with the 
% Pan Tompkins Algorithm

function hrhrvs = hrhrv_pan_tompkin(idx, fs, gr)
    Z = loadEDF(idx);
    K = timetable2table(Z(:,8));

    H = height(K); % the length of the timetable
    C_qrs_amp_raw = cell(H,1);
    C_qrs_i_raw = cell(H,1);
    C_delay = cell(H,1);

    for win = 1:H % window 30 seconds
        Z1 = Z((win:win),:); 
        X = Z1(:,8); % ECG channel (8)
        Y = timetable2table(X);
        Y1 = table2array(Y(:,2));
        Y2 = cell2mat(Y1);
        % Pan tompkins algorithm
        [qrs_amp_raw, qrs_i_raw] = pan_tompkin(Y2,fs,gr);

        C_qrs_amp_raw{win} = qrs_amp_raw;      % qrs_amp_raw per 30 seconds
        C_qrs_i_raw{win} = qrs_i_raw;          % qrs_i_raw per 30 seconds
        % C_delay{win} = delayyy              % delay per 30 seconds
    end

    % Samples to Seconds
    for i = 1:H
        C_qrs_i_raw{i,1} = C_qrs_i_raw{i,1}./fs;
    end

    % Heart rate
    for i = 1:numel(C_qrs_amp_raw)
        heartrateper30sec{i} = numel(C_qrs_amp_raw{i,1});
    end
    heartrateper30sec = heartrateper30sec';
    heartrateper30sec = cell2mat(heartrateper30sec);
    heartrate = heartrateper30sec.*2;
    heartrate = array2table(heartrate);

    % Heart Rate Variability

    for win = 1:H 
        % Determine the RR intervals
        RLocsInterval = diff(C_qrs_i_raw{win,:});
        E{1,win} = num2cell(RLocsInterval); % kathe cell exei ta diasthmata metaksy twn R-R
    end

    F = cell(1,numel(E));

    for b = 1:numel(E)
        F{1,b} = cell2mat(E{1,b});
    end

    for j = 1:numel(F)
        F_std(j) = std(F{1,j}); % Get the std of each column in the matrix.
    end

    HeartRateVariability = F_std;
    HeartRateVariability = HeartRateVariability'; % anastrofos

    HeartRateVariability = array2table(HeartRateVariability);

    hrhrvs = [heartrate HeartRateVariability];
    hrhrvs.Annotations = Z.Annotations;
    
end

