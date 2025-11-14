function [t] = tt_d_markov_s_modified(samples_ind,n, r)
% Compute a rank-r tt-format of an empirical tensor constructed by samples
% Input
%   - samples_ind (given by indices)
%   - n 
%   - r (ranks, scalar or row vector of length d - 1)
% Output
%   - t (tt_tensor object)

d = size(samples_ind,2); % dimension
N = size(samples_ind,1); % sample size
T = cell(d,1); % to store cores 
PT = cell(d-1,1); % to store A

if length(r)==1 % input rank = scalar
    r = r*ones(1,d-1); % [r, r, ..., r]    
end


for k = 1:(d-1)
    if k==1
        % estimate p(x_1, x_2)
        P_hat = zeros(n,n);
        [val,~,ic] = unique(samples_ind(:,[k, k+1]), 'rows');
        freq = groupcounts(ic) / N; % p_hat(x_1, x_2)
        for i = 1:length(freq)
            P_hat(val(i,1),val(i,2)) = freq(i);
        end
        
        if n == r(k)
            B = P_hat;
        else
            % svd p(x1, x2)
            [B,~] = svd(P_hat);
            B = B(:,1:r(k)); % n * r(k)
        end
        
        % update the core
        T{k} = B;
        
        % update A = B
        PT{k} = T{k};     
    else
        % estimate the joint p(x_{k-1], x_k, x_{k+1})
        P_hat = zeros(n,n,n);
        [val,~,ic] = unique(samples_ind(:,[k-1, k, k+1]), 'rows');
        freq = groupcounts(ic) / N; % p_hat(x_{k-1}, x_k, x_{k+1})
        for i = 1:length(freq)
            P_hat(val(i,1),val(i,2),val(i,3)) = freq(i);
        end
        
        P_hat = reshape(P_hat,n^2,n);
        
        if n == r(k)
            B = P_hat;
        else
            % svd p(x_{k-1}, x_{k}; x_{k+1})
            [B,~] = svd(P_hat);
            B = B(:,1:r(k)); % % n^2 * r(k)
        end
        
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
p = reshape(sum(reshape(P_hat,n,n,n),1),n,n);
T{d} = PT{d-1} \ p; % the last core

t = T2t(T,n,d,r); % to a tt_tensor object (TT-Toolbox)
end

