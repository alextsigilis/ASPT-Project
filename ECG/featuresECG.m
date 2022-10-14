% Orfanos Dimitrios, 9579

% Features ECG

function features_ECG = featuresECG(idx, fs, gr)
    Z = loadEDF(idx);
    heartrates = heartrate(idx);
    hrhrvs = hrhrv_pan_tompkin(idx, fs, gr);

    heartrates = heartrates(:,1);
    hrhrvs = hrhrvs(:,2);

    features_ECG = [heartrates hrhrvs];
    features_ECG.Annotations = Z.Annotations;

end
