% Orfanos Dimitrios, 9579

% Heart rate - Heart Rate variability

function new_table = hrhrv(idx,channel,fs)
    Z = loadEDF(idx);
    heartrates = heartrate(idx);
    HeartRateVariability = HRV(idx,channel,fs);
    new_table = [heartrates(:,1) HeartRateVariability(:,1)];
    new_table.Annotations = Z.Annotations;
end





