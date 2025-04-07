beta = 1; lambda = 1; h = 1;    % parameters for GL Potential
a = -4;
b = 4;

% numerical approximation using 
n = 50;
[x_lg, w] = lgwt(n, a, b); % roots of the Lengendre polynomial

x_lg = flip(x_lg); % ascending order
w = flip(w);
W = diag(w); % n * n
W_sq = sqrt(W);

d = 5;
p_lg = gl_cross(x_lg',d,beta,lambda,h,1e-10); % tt_cross of p(x_lg) (n^d)
pw = tmult(W, p_lg); % p(x_lg)*w

% normalize
nconst = sum(pw); % normalizing constant
p_lg = p_lg ./ nconst; % normalize

% L2 norm 
l2_p = norm(tmult(W_sq, p_lg));

MM = [7, 9, 11, 13, 15];
err_a = zeros(1, length(MM));

rep = 20;
err_e = zeros(rep, length(MM));

% basis functions evaluated at x_lg
% fourier
M_max = max(MM);
U_lg = zeros(n, M_max); % (n * M_max)
mid = (a + b) / 2;
hal = (b - a) / 2;
x_lg_std = (x_lg-mid)/hal; % standardize
for j = 1:M_max
    if j==1
        U_lg(:,j) = 1 / sqrt(2*hal);
    elseif rem(j, 2) == 0
        U_lg(:,j) = cos(j/2*pi*x_lg_std) / sqrt(hal);
    else
        U_lg(:,j) = sin((j-1)/2*pi*x_lg_std) / sqrt(hal);
    end
end

r = 3;
n_bin = 300;

rng(12)
N = 1000000;
sigma = 0.5 ; % for MCMC
for i = 1:rep
    disp(['running batch ', num2str(i)])
    
    % sampling
    [samples] = gl_mh(N, d, beta, lambda, h, [a,b], sigma);

    for j = 1:length(MM)
        M = MM(j);
                
        % basis expansion
        U = U_lg(:,1:M); % using M basis functions (n * M)
        p_t = tmult(U'*W, p_lg); % p_tilde(beta) (M^d)
            
        % approximation error
        err_a(j) = norm(tmult(W_sq, p_lg - tmult(U, p_t))) / l2_p;

        % estimation
        p_ts = tt_c_markov_s(samples,M,r,a,b,n_bin,"Fourier");
        err_e(i, j) = norm(p_t - p_ts) / l2_p;
    end
end

save('gl_continuous5.mat')