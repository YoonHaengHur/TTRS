%% n = 5 / Fixed N
supp = [-2, -1, 0, 1, 2];
n = length(supp);
beta = 0.2;

r = 2; 
N = 50000;
dd = 3*(1:10);
rep = 20;
errors1 = zeros(length(dd), rep);
errors2 = zeros(length(dd), rep);

rng(111)
for i = 1:length(dd)
    d = dd(i);
    J = interaction(d);
    disp(['start d = ', num2str(d)])

    p_cross = ising_cross(supp,d,beta,J,1e-14);
    p_cross = p_cross / sum(p_cross); % normalize to a density
    
    for j = 1:rep
        disp(['- running batch ', num2str(j)])
        
        [~,samples_ind] = ising_gibbs(N,supp,d,beta,J);
        [t] = tt_d_markov_s(samples_ind,n,[r,(r+1)*ones(1,d-3),r]);
        errors1(i,j) = norm(t-p_cross) / norm(p_cross);

        [t2] = tt_d_markov_s2(samples_ind,n,[r,(r+1)*ones(1,d-3),r]);
        errors2(i,j) = norm(t2-p_cross) / norm(p_cross);
    end
end

save('ising5_fixed_N.mat')


%%
function [J] = interaction(d)
% J(i, j) = - 1 / (1 + |i-j|) if |i-j|<=2 and 0 otherwise

J = zeros(d, d);
J = J + diag(ones(1,d));
J = J + (1/2) * (diag(ones(1,d-1),-1) + diag(ones(1,d-1),1));
J = J + (1/3) * (diag(ones(1,d-2),-2) + diag(ones(1,d-2),2));
J = -J;
end







