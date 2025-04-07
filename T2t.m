function [t] = T2t(T,n,d,r)
% Input
%   - n (node, scalar)
%   - d (dimension, scalar)
%   - r (ranks, row vector of length d - 1)
% convert T = cell(d, 1) to a tt_tensor object t (TT-Toolbox)
% T{1} = n * r(1)
% T{2} = r(1) * n * r(2)
%       .
%       .
% T{d-1} = r(d-2) * n * r(d-1)
% T{d} = r(d-1) * n
% t = tt_tensor object having T{k} as k-th core 

t = tt_tensor;
t.d = d; % dimension
t.r = [1,r,1]'; % ranks
t.n = n*ones(d,1); % mode sizes 

ps = cumsum([1;t.n.*t.r(1:d).*t.r(2:d+1)]);
t.core = zeros(ps(d+1)-1,1); % vectorize each core and stack all of them

for k=1:d
    cr = T{k};     
    t.core(ps(k):ps(k+1)-1) = cr(:);      
end

t.ps=ps;
end