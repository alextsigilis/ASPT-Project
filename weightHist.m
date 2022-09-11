function weightHist(QPC, f)

% Boolean arrays of sleep stages 
W  = QPC{:, "Annotations"} == "Sleep stage W";
R  = QPC{:, "Annotations"} == "Sleep stage R";
N1 = QPC{:, "Annotations"} == "Sleep stage N1";
N2 = QPC{:, "Annotations"} == "Sleep stage N2";
N3 = QPC{:, "Annotations"} == "Sleep stage N3";

% List of boolean arrays
masks  = {W, R, N1, N2, N3}; 

% List of sleep stages
stages = [              ...
    "Sleep stage W" ,   ... 
    "Sleep stage R" ,   ...
    "Sleep stage N1",   ...
    "Sleep stage N2",   ...
    "Sleep stage N3"]; 

% K: (integer) number of features in QPC table
% M: (integer) number of sleep stages
% N: (integer) number of EEG wave types
K = 9;  M = 5; N = 4;

% F: array of frequency bounds for delta, theta, alpha and beta waves
% coeff: array of normalization factors for weighted histogram
F = [0.0; 4.0; 8.0; 16.0; 32.0];
% coeff = [];

% Estimate normalization coefficients for 
% every region of the bicoherence matrices
% for n = 1:1:4                 
%     for m = n:1:4
%         u = (f > F(n)) & (f <= F(n+1));
%         v = (f > F(m)) & (f <= F(m+1));             
%         w = v & u' & (f >= f') & (f + f' <= F(end));
%         coeff = [coeff sum(w(:))];
%     end
% end

% Discard unnecessary frequencies
if any(f < 0) 
    f = f(f > F(0) & f <= F(end));
end

% Iterate through all sleep stages
for m = 1:1:M                  
    % boolean array of sleep stage Annotations
    idx = masks{m};

    % Histogram array
    h = zeros(numel(f));

    % Iterate through all bicoherence sub-regions
    for k = 1:1:K
        % X: sub-array of QPC features
        X = cell2mat(QPC{idx,k});

        % f1: x-coordinates of histogram
        % f2: y-coordinates of histogram
        % b:  b(f1,f2) "weight" of histogram observations
        f1 = X(:,1);                
        f2 = X(:,2);                
        b  = X(:,3) / sum(X(:,3));  

        % update weighted histogram
        for l = 1:1:size(X,1)
            x = f1(l); y = f2(l); z = b(l);
            h(x,y) = h(x,y) + z;
        end
    end

    % Plot weighted histogram of QPC pairs
    figure(m);
    contourf(f,f,h,16,'LineColor','none');
    
    xlabel("f_1"); 
    ylabel("f_2"); 
    title(sprintf("weighted QPC histograms, %s",stages(m)));

    caxis([0 1]);
    colormap("turbo");
    colorbar;
end

end