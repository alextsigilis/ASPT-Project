% Orfanos Dimitrios, 9579

% ===================================================================
% Estimating and Plotting the Heart Rates and their Histograms
% with the findpeaks
% ===================================================================
%
% Arguments:
% i:  (integer 1-154) ID of the patient 
% channel: (integer 8) the ECG channel
% fs: frequency
%
% ===================================================================
%
% Return Variable: 
% heartrate: array with the Heart Rate per 30 seconds

function [heartrate,HRV] = hr(i,channel,fs)

    % Progress status
    fprintf("Processing patient %d\n\n",i);

    % Make sure the input file exists before attempting to open it
    if ~isfile(sprintf("SN%03d.edf",i))
        fprintf("Patient %d does not exist\n\n",i);
    end

    % Load the ECG recordings and sleep stage labels
    Z = loadEDF(i);
    K = timetable2table(Z(:,8));
    H = height(K); % the length of the timetable
    % hr_30sec: heart rate per 30 seconds
    hr_30sec = zeros(H,1);

    for win = 1:H % window 30 second
        Z1 = Z((win:win),:); 
        X = Z1(:,channel);
        Y = timetable2table(X);
        Y1 = table2array(Y(:,2));
        Y2 = cell2mat(Y1);
        % pks: Local maxima
        % locs: Peak locations
        [pks,locs] = findpeaks(Y2,fs,'MinPeakDistance',0.6,'Annotate','extents');
        hr_30sec(win) = numel(pks);
        
        if locs < 30*win
            % Determine the RR intervals
            RLocsInterval = diff(locs);

            % Derive the HRV signal
            tHRV = locs(2:end);
            HRV = 1./RLocsInterval;
        end
        
%         k = height(pks)
%         l = height(Y2)
%         y = linspace(0,l,k)
%         
%         % Plotting the peaks
%         figure
%         plot(Y2, '-b')
%         hold on
%         plot(y, pks, 'pg', 'MarkerFaceColor','g')
%         hold off
%         grid
    end
        
    % heartrate: heart rate per 1 minute
    heartrate = hr_30sec*2;
    
    % Plots
    figure;
    plot(heartrate)
    xlabel('Time')
    title('Heart rate')
    
    figure;
    hist(heartrate,15)
    title('Histogram of the Heart rate')

end

