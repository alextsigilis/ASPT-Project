% Orfanos Dimitrios, 9579

% Histograms for the heart rate

function [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5] = histogramhr(heartrate, nbins)
    s1 = heartrate.Annotations == "Sleep stage W";
    [y1, x1] = hist(heartrate.heartrate(s1), nbins);
    y1 = y1 / sum(y1);

    s1 = heartrate.Annotations == "Sleep stage N1";
    [y2, x2] = hist(heartrate.heartrate(s1), nbins);
    y2 = y2 / sum(y2);

    s1 = heartrate.Annotations == "Sleep stage N2";
    [y3, x3] = hist(heartrate.heartrate(s1), nbins);
    y3 = y3 / sum(y3);

    s1 = heartrate.Annotations == "Sleep stage N3";
    [y4, x4] = hist(heartrate.heartrate(s1), nbins);
    y4 = y4 / sum(y4);

    s1 = heartrate.Annotations == "Sleep stage R";
    [y5, x5] = hist(heartrate.heartrate(s1), nbins);
    y5 = y5 / sum(y5);

end