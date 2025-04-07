function [samples,samples_ind] = gl_gibbs(N, supp, d, beta, lambda, h)
% Gibbs sampler from the discrete GL
% Input
%   - N (# of samples)
%   - supp (support of each coordinate = n values)
%   - d (dimension)
%   - beta, lambda, h (parameters for GL potential) 
% Output 
%   - samples (exact values)
%   - samples_ind (by indices)

n = length(supp);

% draw 2N samples and discard the first half (burn-in)
samples = zeros(2*N,d); % each row represents a sample
samples_ind = zeros(2*N, d); 

% samples_ind represents samples via indices:
% ex) (1, 2, 3) in samples_ind <-> (supp(1), supp(2), supp(3)) in samples


X0 = zeros(1,d); % initial sample
samples(1,:) = X0; % set the first row which will be updated in the first iteration
for i = 1:(2*N)
    for j = 1:d
        if j == 1
            prob = exp(-beta*(lambda/2*supp.^2/h^2 + lambda/2*(supp-samples(i,j+1)).^2/h^2 + 1/(4*lambda)*(supp.^2-1).^2)); % conditional probabilities of j-th coord = supp(1)~supp(n)
        elseif j == d
            prob = exp(-beta*(lambda/2*(samples(i,j-1)-supp).^2/h^2 + lambda/2*supp.^2/h^2 + 1/(4*lambda)*(supp.^2-1).^2));
        else            
            prob = exp(-beta*(lambda/2*(samples(i,j-1)-supp).^2/h^2 + lambda/2*(supp-samples(i,j+1)).^2/h^2 + 1/(4*lambda)*(supp.^2-1).^2));
        end
        
        prob = prob / sum(prob); % normalize to a prob vector    
        samples_ind(i, j) = randsrc(1, 1, [1:n; prob]); % draw index using the prob vector (randsrc from communitcations toolbox)
        samples(i,j) = supp(samples_ind(i, j)); % update the j-th coordinate with the exact value from the support
    end
    if i < 2*N
        samples(i+1,:) = samples(i,:); % set the next row with the current sample which will be updated in the following iteration
    end
end

samples = samples((N+1):end,:); % burn-in
samples_ind = samples_ind((N+1):end,:);


end

