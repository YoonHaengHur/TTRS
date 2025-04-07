function [t_new] = tmult(X, t)
% X = m * n
% t = tt_object with k-th core of size t.r(k) * t.n(k) (= n) * t.r(k+1)
% contract each core G_k and X so that the size be t.r(k) * m * t.r(k+1)

m = size(X, 1);
d = t.d;

t_new = tt_tensor;
t_new.d = d;
t_new.n = m*ones(d,1); % mode sizes 
t_new.r = t.r;

ps = cumsum([1;t_new.n.*t_new.r(1:d).*t_new.r(2:d+1)]);
t_new.core = zeros(ps(d+1)-1,1); % vectorize each core and stack all of them

for k=1:d
    cr = t.core(t.ps(k):t.ps(k+1)-1);
    cr = reshape(cr,t.r(k),t.n(k),t.r(k+1)); % r1 * n * r2
    cr = cmult(X,cr); % r1 * m * r2
    t_new.core(ps(k):ps(k+1)-1) = cr(:);      
end

t_new.ps=ps;
end


 
function [Z] = cmult(X, Y)
% X = m * n 
% Y = r1 * n * r2
% Z = contraction of X and Y such that the size = r1 * m1 * r2

m = size(X, 1);
n = size(X, 2);
r1 = size(Y, 1);
r2 = size(Y, 3);
tmp = X*reshape(permute(Y,[2,1,3]),n,[]); % m * (r1 * r2)
Z = permute(reshape(tmp, m, r1, r2),[2,1,3]);
end

