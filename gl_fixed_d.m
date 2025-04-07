%% Fixed d
% support of each coordinate (n values for each coordinate)
n = 9;
supp = linspace(-4,4,n);
d = 8;
beta = 1; lambda = 1; h = 1; % parameters for GL Potential
p_cross = gl_cross(supp,d,beta,lambda,h,1e-14); 
p_cross = p_cross / sum(p_cross); % normalize to a density

r = 3; % TT-ranks

NN = 2.^(8:17);
rep = 20;
errors = zeros(length(NN), rep);

rng(111)
for i = 1:length(NN)
    N = NN(i);
    disp(['start N = ', num2str(N)])

    for j = 1:rep
        disp(['- running batch ', num2str(j)])

        [~,samples_ind] = gl_gibbs(N,supp,d,beta,lambda,h);
        [t] = tt_d_markov_s(samples_ind,n,r); 
        errors(i, j) = norm(t-p_cross) / norm(p_cross);
    end
end

save('gl_fixed_d.mat')

