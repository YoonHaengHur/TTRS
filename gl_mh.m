function [samples] = gl_mh(N, d, beta, lambda, h, supp, sigma)
% Metropolis-Hastings algorithm to sample from GL potential 
% Input
%   - N (# of samples)
%   - d (dimension)
%   - beta, lambda, h (parameters for GL potential)
%   - supp = [a, b] so that each coordinate is contained in [a, b]
%   - sigma (variance of the proposal increment distribution) 
% Output 
%   - samples 

% draw 2N samples and discard the first half (burn-in)
samples = zeros(2*N,d); % each row represents a sample 

% proposal density (Markov kernel)
p_proposal = @(x,y) mvnpdf(y, x, sigma^2*eye(d));

x = (supp(2) - supp(1))*rand(1, d) + supp(1);
for i=1:(2*N)
    x_new = x + sigma*mvnrnd(zeros(1,d),eye(d));
    alpha_1 = p(x_new, beta, lambda, h, supp) * p_proposal(x_new,x); 
    if alpha_1 == 0
        alpha = 0;
    else
        alpha_2 = p(x, beta, lambda, h, supp) * p_proposal(x,x_new);
        alpha = alpha_1 / alpha_2;
    end
    if rand <= min(alpha,1) 
        samples(i,:) = x_new; % accept
    else 
        samples(i,:) = x; % reject   
    end
    x=samples(i,:);
end

samples = samples((N+1):end,:); % burn-in

end

% target density (x is a row vector)
function [p] = p(x, beta, lambda, h, supp)
    if (supp(1)<=x)&(x<=supp(2))
        p = exp(-beta*(lambda/2*sum(([0,x]-[x,0]).^2)/h^2+1/(4*lambda)*sum(([0,x].^2-1).^2))); 
    else 
        p = 0;
    end
end
