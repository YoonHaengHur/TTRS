function [p_cross] = ising_cross(supp, d, beta, J, epsilon)

n = length(supp);
% density interms of indices
p_ind = @(x) exp(-beta*(supp(x)*J*supp(x)')); 
   
p_cross = dmrg_cross(d,n,p_ind,epsilon); % cross
end

