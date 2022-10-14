% Orfanos Dimitrios, 9579

% Function Description:
% This function exports and plots the heart rate.
% ---------------------------------------------------------
%
% Function Arguments:
% idx: (integer) the index of the patient's heart rate 
%
% Return Variables:
% heartrate: (array) An array with the values of Heart Rate
% per 30 seconds and the "Annotations" column which contains 
% the sleep stage labels for every segment.
% =========================================================

function heartrate = heartrate(idx)
    X = loadEDF_HeartRateECG(idx);
    Y = timetable2table(X);
    Y1 = table2array(Y(:,2));
    Y2 = cell2mat(Y1);
    
    % Plotting the continuous heart rate
    figure;
    plot(Y2)
    title("Continuous Heart Rate")

    for i = 1:numel(Y1)
        heartrate(i) = mean(Y1{i,1}); % Get the mean of each column in the matrix.
    end
    heartrate = heartrate';
    
    % Plotting the heart rate per 30 seconds
    figure;
    plot(heartrate)
    title("Heart Rate per 30 seconds")
    
    heartrate = array2table(heartrate);
    heartrate.Annotations = X.Annotations;
    
    [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5] = histogramhr(heartrate, 20);
    % Display the histograms
    figure; hold on; grid on;
    plot(x1,y1,'r',x2,y2,'b',x3,y3,'k',x4,y4,'g',x5,y5,'y');
    xlabel("Heart rate");
    ylabel("Probability");
    title("Histogram of the Heart Rate");
    legend("Wake", "N1", "N2", "N3", "REM");
end