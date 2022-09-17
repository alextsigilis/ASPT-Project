% ========================================================================
% QPC features:
% 1) Number of bicoherence peaks
% 2) Average bicoherence values of peaks
% 3) Average euclidean distance of the peaks from the axis
% 4) Euclidian distance of the most distant peak from the axis origin
% 5) Center of mass
% ========================================================================

function Y = QPCfeatures(QPC, f)
    N = size(QPC,1);
    
    names = [           ...
        "numOfPeaks",   ...
        "avgBic",       ...
        "avgDist",      ...
        "maxDist",      ...
        "center",       ...
        "Annotations"];

    types = [repmat("double",1,numel(names)-1) "string"];
    sz = [N numel(names)];

    Y = table(                      ...
        'Size',          sz,        ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    for i = 1:1:N
        % array of bicoherence peaks
        peaks = cell2mat(QPC{i,1});

        % locations and values of bicoherence peaks
        x = f(peaks(:,1)); 
        y = f(peaks(:,2));
        b = peaks(:,3); b = b(:);

        if numel(x) == 0
            Y{i,"numOfPeaks"} = nan;
            Y{i,"avgBic"} = nan;
            Y{i,"avgDist"} = nan;
            Y{i,"maxDist"} = nan;
            Y{i,"center"} = nan;
            continue;
        end

        % Normalize x and y
        x = x(:) / f(end);
        y = y(:) / f(end);

        % distances of the bicoherence peaks from the origin
        % d1 = x + y;
        d2 = sqrt(x.^2 + y.^2);

        % number of bicoherence peaks
        Y{i,"numOfPeaks"} = numel(b(b>=0.01));

        % average of peak bicoherence
        Y{i,"avgBic"} = sum(d2 .* b.^(1/4)) / sum(d2);

        % (weighted) average of euclidean distance of bicoherence
        % peaks from the origin of the frequency axes
        Y{i,"avgDist"} = sum(d2 .* b.^(1/4));

        % (weighted) maximum euclidean distance of bicoherence
        % peaks from the origin of the frequency axes.
        Y{i,"maxDist"} = max(d2 .* b.^(1/4));

        % center of mass
        r = 1.2;
        cx = sum(b.^(1/2) .* x) / sum(b.^(1/2));
        cy = sum(b.^(1/2) .* y) / sum(b.^(1/2));
        Y{i,"center"} = (cx^r + cy^r)^(1/r);
    end

    Y.Annotations = QPC.Annotations;
end