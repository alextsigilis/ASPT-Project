% ========================================================================
% QPC features:
% 1) Number of bicoherence peaks
% 2) Average bicoherence values of peaks
% 3) Average euclidean distance of the peaks from the axis
% 4) Euclidian distance of the most distant peak from the axis origin
% 5) Center of mass
% ========================================================================

function Y = QPCfeatures(QPC, f, weighted)
    N = size(QPC,1);
    
    names = [           ...
        "numOfPeaks",   ...
        "sumBic",       ...
        "avgDist",      ...
        "maxDist",      ...
        "center",       ...
        "Annotations"];

    types = [           ...
        "double",       ...
        "double",       ...
        "double",       ...
        "double",       ...
        "double",       ...
        "string"];

    sz = [N numel(names)];

    Y = table(                      ...
        'Size',          sz,        ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    if weighted == false
        for i = 1:1:N
            % array of bicoherence peaks
            peaks = cell2mat(QPC{i,1});

            if size(peaks,1) == 0
                Y{i,"numOfPeaks"} = nan;
                Y{i,"sumBic"} = nan;
                Y{i,"avgDist"} = nan;
                Y{i,"maxDist"} = nan;
                Y{i,"center"} = nan;
                continue;
            end
            
            % locations and values of bicoherence peaks
            x = f(peaks(:,1)); x = x(:) / f(end);
            y = f(peaks(:,2)); y = y(:) / f(end);
            b = peaks(:,3);    b = b(:); b = sqrt(b);

            % euclidean distances of bicoherence peaks from the origin
            d = (x.^2 + y.^2).^(1/2); d = d(:);

            % number of bicoherence peaks
            Y{i,"numOfPeaks"} = numel(b);

            % cumulative sum of peak-bicoherence 
            Y{i,"sumBic"} = sum(b);

            % average euclidean distance of bicoherence peaks
            % from the origin of the frequency axes.
            Y{i,"avgDist"} = mean(d);

            % maximum euclidean distance of bicoherence peaks
            % from the origin of the frequency axes.
            Y{i,"maxDist"} = max(d);

            % center of mass
            cx = sum(b .* x) / sum(b); 
            cy = sum(b .* y) / sum(b);
            Y{i,"center"} = sqrt(cx ^ 2 + cy ^ 2);
        end
    
    elseif weighted == true
        for i = 1:1:N
            % array of bicoherence peaks
            peaks = cell2mat(QPC{i,1});

            if size(peaks,1) == 0
                Y{i,"numOfPeaks"} = nan;
                Y{i,"sumBic"} = nan;
                Y{i,"avgDist"} = nan;
                Y{i,"maxDist"} = nan;
                Y{i,"center"} = nan;
                continue;
            end

            % locations and values of bicoherence peaks
            x = f(peaks(:,1)); x = x(:) / f(end);
            y = f(peaks(:,2)); y = y(:) / f(end);
            b = peaks(:,3);    b = b(:); b = sqrt(b);

            % euclidean distances of the bicoherence peaks from the origin
            d = (x.^2 + y.^2).^(1/2); d = d(:);

            % number of bicoherence peaks
            Y{i,"numOfPeaks"} = numel(b);

            % average of peak bicoherence
            Y{i,"sumBic"} = mean(b);

            % (weighted) average of euclidean distance of bicoherence
            % peaks from the origin of the frequency axes
            Y{i,"avgDist"} = sum(d .* b);

            % (weighted) maximum euclidean distance of bicoherence 
            % peaks from the origin of the frequency axes.
            Y{i,"maxDist"} = max(d .* b);
            
            % center of mass
            cx = sum(b .* x) / sum(b); 
            cy = sum(b .* y) / sum(b);
            Y{i,"center"} = sqrt(cx ^ 2 + cy ^ 2);
        end
    end

    Y.Annotations = QPC.Annotations;
end