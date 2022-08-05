function [delta alpha beta theta] = mraEEG(X,idx)
    % fs: (integer) sampling rate of X{:,idx} in Hertz
	% w:  (integer) duration of labeled segments in seconds
	% K:  (integer) number of labeled segments
	% N:  (integer) number of samples per labeled segment
	fs = 256; w = 30; N = w * fs; [K,~] = size(X); 
	 
	x = X{:,idx};
	x = cell2mat(x);
	x = reshape(x, [N K]);		
	
	wt  = modwt(x,"db2",5);
	mra = modwtmra(wt,"db2");
		
	delta = nan(K,3);
	alpha = nan(K,3);
	beta  = nan(K,3);
	theta = nan(K,3);
	
	delta(:,1) = std(mra(6,:,:),0,2);
	delta(:,2) = skewness(mra(6,:,:),0,2);
	delta(:,3) = kurtosis(mra(6,:,:),0,2);
	
	alpha(:,1) = std(mra(5,:,:),0,2);
	alpha(:,2) = skewness(mra(5,:,:),0,2);
	alpha(:,3) = kurtosis(mra(5,:,:),0,2);
	
	beta(:,1)  = std(mra(4,:,:),0,2);
	beta(:,2)  = skewness(mra(4,:,:),0,2);
	beta(:,3)  = kurtosis(mra(4,:,:),0,2);
	
	theta(:,1) = std(mra(3,:,:),0,2);
	theta(:,2) = skewness(mra(3,:,:),0,2);
	theta(:,3) = kurtosis(mra(3,:,:),0,2);
	
	names = {'var' 'skw' 'krt'};
	
	delta = array2table(delta, 'VariableNames', names);
	alpha = array2table(alpha, 'VariableNames', names);
	beta  = array2table(beta,  'VariableNames', names);
	theta = array2table(theta, 'VariableNames', names);
	
	delta = addvars(delta, X.Annotations, 'NewVariableNames', 'Annotations');
	alpha = addvars(alpha, X.Annotations, 'NewVariableNames', 'Annotations');
	beta  = addvars(beta,  X.Annotations, 'NewVariableNames', 'Annotations');
	theta = addvars(theta, X.Annotations, 'NewVariableNames', 'Annotations');
end