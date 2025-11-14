function [t] = tt_marginal_s2(samples_ind,n)
% Compute a rank-n tt-format of an empirical tensor constructed by samples
% Input
%   - samples_ind (given by indices)
%   - n 
% Output
%   - t (tt_tensor object)

d = size(samples_ind,2); % dimension
N = size(samples_ind,1); % sample size
T = cell(d,1); % to store cores G

for k = 1:(d-1)
    % 1. estimate the joint p(x_k, x_{k+1}) and marginal p(x_k)
    P_hat = zeros(n,n);
    p_hat = zeros(n,1);
    [val,~,ic] = unique(samples_ind(:,[k, k+1]), 'rows');
    freq = groupcounts(ic) / N; % p_hat(x_k, x_{k+1})
    for i = 1:length(freq)
        P_hat(val(i,1),val(i,2)) = freq(i);
        p_hat(val(i,1)) = p_hat(val(i,1)) + freq(i);            
    end

    if k==1
        % p(x_1, x_2)
        T{k} = P_hat;            
    else
        % 2. take the inverse of p(x_k)
        p_hat_inv = zeros(n,1);
        pos = p_hat>0;
        p_hat_inv(pos) = 1./p_hat(pos);

        % 3. estimate p(x_{k+1} | x_k)
        P_hat = p_hat_inv .* P_hat;

        % 4. add a leg G_k(arg1, arg2, arg3) with arg1 = arg2
        tmp = zeros(n,n,n);
        for i=1:n
            tmp(i,i,:) = P_hat(i,:); % add a leg by identity
        end        
        T{k} = tmp;
    end
end

T{d} = eye(n);

t = T2t(T,n,d,n*ones(1,d-1)); % to a tt_tensor object (TT-Toolbox)
end

