%% Fixed N
n = 9;
supp = linspace(-4,4,n);
beta = 1; lambda = 1; h = 1; % parameters for GL Potential

N = 50000;
r = 3; % TT-ranks

dd = 3*(1:10);
rep = 20;
errors = zeros(length(dd), rep);

rng(23)
for i = 1:length(dd)
    d = dd(i);
    disp(['start d = ', num2str(d)])

    p_cross = gl_cross(supp,d,beta,lambda,h,1e-10); 
    p_cross = p_cross / sum(p_cross); % normalize to a density
    
    for j = 1:rep
        disp(['- running batch ', num2str(j)])
        
        [~,samples_ind] = gl_gibbs(N,supp,d,beta,lambda,h);    
        [t] = tt_d_markov_s(samples_ind,n,r); 
        errors(i, j) = norm(t-p_cross) / norm(p_cross);
    end
end

save('gl_fixed_N.mat')

