function [t] = tt_marginal_s(samples_ind,n)
% Compute a rank-n tt-format of an empirical tensor constructed by samples
% Input
%   - samples_ind (given by indices)
%   - n 
% Output
%   - t (tt_tensor object)

d = size(samples_ind,2); % dimension
N = size(samples_ind,1); % sample size
T = cell(d,1); % to store cores G

for k = 1:d
    if k==1
        % estimate the marginal p(x_1)
        p1_hat = zeros(n,1);
        x = samples_ind(:,k);
        [val,~,ic] = unique(x, 'sorted');
        freq = groupcounts(ic) / N;
        
        for i=1:length(freq)
            p1_hat(val(i)) = freq(i);
        end
        
        T{k} = diag(p1_hat);
    else
        % 1. estimate the joint p(x_{k-1}, x_k) and marginal p(x_{k-1})
        P_hat = zeros(n,n);
        p_hat = zeros(n,1);
        [val,~,ic] = unique(samples_ind(:,[k-1, k]), 'rows');
        freq = groupcounts(ic) / N; % p_hat(x_{k-1}, x_k)
        for i = 1:length(freq)
            P_hat(val(i,1),val(i,2)) = freq(i);
            p_hat(val(i,1)) = p_hat(val(i,1)) + freq(i);            
        end
        
        % 2. take the inverse of p(x_{k-1})
        p_hat_inv = zeros(n,1);
        pos = p_hat>0;
        p_hat_inv(pos) = 1./p_hat(pos);
        
        % 3. estimate p(x_k | x_{k-1})
        P_hat = p_hat_inv .* P_hat; % p_hat(x_k | x_{k-1})
        
        if k==d
            T{k} = P_hat;
        else
            tmp = zeros(n,n,n);
            for i=1:n
                tmp(:,i,i) = P_hat(:,i); % add a leg by identity
            end        
            T{k} = tmp;
        end
    end    
end

t = T2t(T,n,d,n*ones(1,d-1)); % to a tt_tensor object (TT-Toolbox)
end

