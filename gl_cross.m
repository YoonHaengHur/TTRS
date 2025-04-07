function [p_cross] = gl_cross(supp, d, beta, lambda, h, epsilon)

n = length(supp);

% density interms of indices
p_ind = @(x) exp(-beta*(lambda/2*sum(([0,supp(x)]-[supp(x),0]).^2)/h^2+1/(4*lambda)*sum(([0,supp(x)].^2-1).^2))); 
   
p_cross = dmrg_cross(d,n,p_ind,epsilon); % cross
end

