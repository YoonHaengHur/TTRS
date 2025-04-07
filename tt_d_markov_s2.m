function [t] = tt_d_markov_s2(samples_ind,n,r)
% Compute a rank-r tt-format of an empirical tensor constructed by samples
% Input
%   - samples_ind (given by indices)
%   - n 
%   - r (ranks, scalar or row vector of length d - 1)
% Output
%   - t (tt_tensor object)

d = size(samples_ind,2); % dimension
T = cell(d,1); % to store cores G
PT = cell(d-1,1); % to store A

if length(r)==1 % input rank = scalar
    r = r*ones(1,d-1); % [r, r, ..., r]    
end

% group by unique samples and compute frequencies (empirical tensor)
[s,~,is] = unique(samples_ind,'rows'); % a consists of unique samples
N = size(s,1); % # of unique samples
P = histcounts(is,N,'Normalization','probability'); % frequencies
P = P(:); % vectorized empirical tensor

for k=1:(d-1)
    if k==1
        % p(x1, x2, x3)
        [z,~,iz] = unique(s(:,[k+1,k+2]),'rows');
        
        Y = zeros(n,size(z,1)); % n * n2 where n2 <= n^2
        for y = 1:n
            ind = find(s(:,k) == y);
            l = length(ind); % l could be 0 (zero row)
            if l>0
                for j = 1:l
                    Y(y,iz(ind(j))) = Y(y,iz(ind(j))) + P(ind(j));
                end
            end
        end
        
        % svd p(x1, x2)
        [B,~] = svd(Y);
        B = B(:,1:r(k)); % n * r(k)
        
        % update the core
        T{k} = B;
        
        % update A = B
        PT{k} = T{k};
    elseif k<(d-1)
        % p(x_{k-1}; x_{k}, x_{k+1}, x_{k+2})
        [z,~,iz] = unique(s(:,[k,k+1,k+2]),'rows');
        
        Y = zeros(n,size(z,1)); % n * n2 where n2 <= n^3
        for y = 1:n
            ind = find(s(:,k-1) == y);
            l = length(ind);
            if l>0
                for j = 1:l
                    Y(y,iz(ind(j))) = Y(y,iz(ind(j))) + P(ind(j));
                end
            end
        end
        
        % complete the matrix (n * n^3) before reshaping it
        p = zeros(n,n^3);
        cid = z(:,1) + n*(z(:,2)-1) + n^2*(z(:,3)-1);
        
        % column ordering (x_{k}, x_{k+1}, x_{k+2})
        % (1, 1, 1) = 1
        % (2, 1, 1) = 2    
        %   .
        %   .
        %   .
        % (n, 1, 1) = n
        % (1, 2, 1) = n+1
        % (2, 2, 1) = n+2
        %   .
        %   .
        %   .
        % (n, n, 1) = n^2
        % (1, 1, 2) = n^2 + 1
        
        p(:,cid) = Y;
        
        % svd p(x_{k-1}, x_{k}; x_{k+1}, x_{k+2})
        [B,~] = svd(reshape(p,n^2,n^2));
        B = B(:,1:r(k)); % n^2 * r(k)
        
        % solve to obtain a core
        B = reshape(B,n,n*r(k)); % n * nr(k)
        G =  PT{k-1} \ B; % r(k-1) * nr(k)
        
        % update the core
        T{k} = reshape(G,r(k-1),n,r(k)); 
        
        % update A = sum_{x_{k-1}} B(x_{k-1}; x_{k}, alpha_{k})
        PT{k} = reshape(sum(B),n,r(k));
    else % k == d-1
        % p(x_{k-1}; x_{k}, x_{k+1})
        [z,~,iz] = unique(s(:,[k,k+1]),'rows');
        
        Y = zeros(n,size(z,1)); % n * n2 where n2 <= n^2
        for y = 1:n
            ind = find(s(:,k-1) == y);
            l = length(ind);
            if l>0
                for j = 1:l
                    Y(y,iz(ind(j))) = Y(y,iz(ind(j))) + P(ind(j));
                end
            end
        end
        
        % complete the matrix (n * n^2) before reshaping it
        p = zeros(n,n^2);
        cid = z(:,1) + n*(z(:,2)-1);
        
        % column ordering (x_{k}, x_{k+1})
        % (1, 1) = 1
        % (2, 1) = 2    
        %   .
        %   .
        %   .
        % (n, 1) = n
        % (1, 2) = n+1
        % (2, 2) = n+2
        %   .
        %   .
        %   .
        % (n, n) = n^2
        
        p(:,cid) = Y;
        
        % svd p(x_{k-1}, x_{k}; x_{k+1})
        [B,~] = svd(reshape(p,n^2,n));
        B = B(:,1:r(k)); % n^2 * r(k)
        
        % solve to obtain a core
        B = reshape(B,n,n*r(k)); % n * nr(k)
        G =  PT{k-1} \ B; % r(k-1) * nr(k)
        
        % update the core
        T{k} = reshape(G,r(k-1),n,r(k)); 
        
        % update A = sum_{x_{k-1}} B(x_{k-1}; x_{k}, alpha_{k})
        PT{k} = reshape(sum(B),n,r(k));
    end
end

% p(x_{d-1}, x_{d})
p = reshape(sum(p),n,n);
T{d} = PT{d-1} \ p; % the last core

t = T2t(T,n,d,r); % to a tt_tensor object (TT-Toolbox)
end

