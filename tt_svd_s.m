function [t] = tt_svd_s(samples_ind,n,r)
% Compute a rank-r tt-format of an empirical tensor constructed by samples
% Input
%   - samples_ind (given by indices)
%   - n 
%   - r (ranks, scalar or row vector of length d - 1)
% Output
%   - t (tt_tensor object)

d = size(samples_ind,2); % dimension
T = cell(d,1); % to store cores G

% group by unique samples and compute frequencies
[a, ~, b] = unique(samples_ind, 'rows'); % a consists of unique samples
N = size(a, 1); % # of unique samples
x = histcounts(b, N, 'Normalization', 'probability'); % frequencies
weights = x'; % convert to a column vector aka probability vector 

%{
% the previous step may be skipped
% but can be slow if there are many overlapping samples
a = samples_ind;
N = size(a, 1);
weights = ones(N, 1) / N;
%}

for i=1:(d-1)
    
    % group by x_{i+1} ~ x_{d}
    [tmp1,~,tmp2] = unique(a(:,(i+1):end), 'rows');
    
    % reshape weights (alpha_{i-1}, x_i) by (x_{i+1}, ..., x_{d})
    if i==1
        Y = zeros(n, size(tmp1,1));
        for k = 1:n
            ind = find(a(:,i) == k);
            m = length(ind);
            if m>0
                for j = 1:m
                    Y(k,tmp2(ind(j))) = Y(k,tmp2(ind(j))) + weights(ind(j),:);
                end
            end
        end
    else
        Y = zeros(r*n, size(tmp1,1));
        for k = 1:n
            ind = find(a(:,i) == k);
            m = length(ind);
            if m>0
                for j=1:m
                    Y((k-1)*r+(1:r),tmp2(ind(j))) = Y((k-1)*r+(1:r),tmp2(ind(j))) + weights(ind(j),:)'; 
                end
            end
        end
    end
    
    % Finding r left singular vectors of Y
    % Eigendecompose YY' instead of SVD(Y) (much faster)
        
    % eigen decomposition
    A = Y * Y';
    [U,~] = eigs(A, r);
    
    % update the core
    if i==1
        T{i} = reshape(U, n, r); 
    else
        T{i} = reshape(U, r, n, r); 
    end
    
    % update the weights 
    weights_new = zeros(N, r);
    if i==1
        for k = 1:n
            ind = (a(:,i) == k);
            m = sum(ind);

            if m>0
                weights_new(ind,:) = weights(ind,:) * U(k,:); 
            end
        end       
    else
        for k = 1:n
            ind = (a(:,i) == k);
            m = sum(ind);

            if m>0
                weights_new(ind,:) = weights(ind,:) * U((k-1)*r+(1:r),:); 
            end
        end
    end
    weights = weights_new; 
end
lastcore = zeros(r, n);

for k = 1:n
    ind = (a(:,d) == k);
    m = sum(ind);
    
    if m>0
        lastcore(:,k) = weights(ind,:)' * ones(m,1); 
    end
end
T{d} = lastcore; % the last core


% to a tt_tensor object (TT-Toolbox)
t = tt_tensor;
t.d  = d; % dimension
t.r  = [1, r*ones(1,d-1), 1]'; % ranks
t.n = n*ones(d,1); % mode sizes 

ps=cumsum([1;t.n.*t.r(1:d).*t.r(2:d+1)]);
t.core=zeros(ps(d+1)-1,1); % vectorize each core and stack all of them
t.core(ps(1):ps(2)-1)=T{1}(:);

for i=1:d
    cr=T{i};     
    t.core(ps(i):ps(i+1)-1) = cr(:);      
end

t.ps=ps;

end
