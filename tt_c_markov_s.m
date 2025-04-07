function [t] = tt_c_markov_s(samples,M,r,a,b,n,basis)
% Compute a rank-r tt-format of an empirical tensor constructed by samples
% Input
%   - samples
%   - M (# of basis functions = Legendre polynomials)
%   - r (ranks, scalar or row vector of length d - 1)
%   - [a, b] 
%   - n (for binning the samples)
%   - basis (Fourier or Legendre)
% Output
%   - t (tt_tensor object)

m = (a + b) / 2;
h = (b - a) / 2;
c = 1 / sqrt(2 * h); % phi_{1}(x) = sqrt(1 / (b-a)) = 1 / sqrt(2h)
d = size(samples,2); % dimension
T = cell(d,1); % to store cores G
PT = cell(d-1,1); % to store A

if length(r)==1 % input rank = scalar
    r = r*ones(1,d-1); % [r, r, ..., r]    
end

% group by unique samples and compute frequencies (empirical tensor)
[s,~,is] = unique(samples,'rows'); % a consists of unique samples
N = size(s,1); % # of unique samples
P = histcounts(is,N,'Normalization','probability'); % frequencies
P = P(:); % vectorized empirical tensor

x = linspace(a,b,n); % grid points to bin the samples
U = zeros(n, M);
x_std = (x - m) / h;
if basis == "Fourier"
    for beta = 1:M
        if beta==1
            U(:,beta) = 1 / sqrt(2*h);
        elseif rem(beta, 2) == 0
            U(:,beta) = cos(beta/2*pi*x_std) / sqrt(h);
        else
            U(:,beta) = sin((beta-1)/2*pi*x_std) / sqrt(h);
        end
    end
elseif basis == "Legendre"
    for beta = 1:M
        U(:,beta) = legendreP(beta-1,x_std) * sqrt((2*beta-1)/2) / sqrt(h);
    end
end

for k=1:d
    if k==1
        % p(x1, x2)
        [y,~,iy] = unique(s(:,[k,k+1]),'rows');
        ny = size(y,1);
        Y = zeros(ny,1); % ny * 1
        for i = 1:N
            Y(iy(i)) = Y(iy(i)) + P(i);
        end
        
        bY = bin2(y,x); % ny * 2
        
        Y2 = zeros(n,n);
        for i = 1:ny
            i1=bY(i,1);
            i2=bY(i,2);
            Y2(i1,i2) = Y2(i1,i2) + Y(i);
        end
                
        % svd p_t(beta1, gamma1)
        Y = U' * Y2 * U ; % M * M
        Y = reshape(Y,M,M);
        [B,~] = svd(Y);
        B = B(:,1:r(k)); % M * r(k)
        
        % update the core
        T{k} = B;
        
        % update A = B
        PT{k} = T{k};
    elseif k<d    
        % p(x_{k-1}; x_{k}, x_{k+1})
        [y,~,iy] = unique(s(:,[k-1,k,k+1]),'rows');
        ny = size(y,1);
        Y = zeros(ny,1); 
        for i = 1:N
            Y(iy(i)) = Y(iy(i)) + P(i);
        end
        
        bY = bin3(y,x); % ny * 3
        
        Y2 = zeros(n,n,n);
        for i = 1:ny
            i1=bY(i,1);
            i2=bY(i,2);
            i3=bY(i,3);
            Y2(i1,i2,i3) = Y2(i1,i2,i3) + Y(i);
        end
        
        Y2 = reshape(Y2, n, n^2);
        
        % p_t(beta_{k-1},j_{k};gamma_{k})
        Y = U' * Y2 * kron(U, U); % M * M^2
        Y = reshape(Y, M^2, M); % M^2 * M
        
        if k>2
            Y = c^(k-2) .* Y;
        end
        
        [B,~] = svd(Y);
        B = B(:,1:r(k)); % M^2 * r(k)
        
        % solve to obtain a core
        B = reshape(B,M,M*r(k)); % M * Mr(k)
        G =  PT{k-1} \ B; % r(k-1) * Mr(k)
        
        % update the core
        T{k} = reshape(G,r(k-1),M,r(k)); 
        
        % update A = B(1; j_{k};gamma_{k})
        PT{k} = reshape(B(1,:),M,r(k)); 
    else
        % p(x_{d-1}, x_{d})        
        [y,~,iy] = unique(s(:,[k-1,k]),'rows');
        ny = size(y,1);
        Y = zeros(ny,1); % ny * 1
        for i = 1:N
            Y(iy(i)) = Y(iy(i)) + P(i);
        end
        
        bY = bin2(y,x); % ny * 2
        
        Y2 = zeros(n,n);
        for i = 1:ny
            i1=bY(i,1);
            i2=bY(i,2);
            Y2(i1,i2) = Y2(i1,i2) + Y(i);
        end
        
        % p_t(beta_{d-1},j_{d})
        Y = U' * Y2 * U; % M * M
        
        if k>2
            Y = c^(k-2) .* Y;
        end
        
        % update the core
        T{k} = PT{k-1} \ Y; % r * M        
    end
end

t = T2t(T,M,d,r); % to a tt_tensor object (TT-Toolbox)
end


function [id] = bin3(s, x)
% x consists of n points
% s = N * 3

s1=s(:,1);
s2=s(:,2);
s3=s(:,3);
id=[bin1(s1,x), bin1(s2,x), bin1(s3,x)];
end

function [id] = bin2(s, x)
% x consists of n points
% s = N * 2

s1=s(:,1);
s2=s(:,2);
id=[bin1(s1,x), bin1(s2,x)];
end

function [id] = bin1(s,x)
% x consists of n points
% s = N * 1
% id = N * 1
    [~,id] = min(abs(s(:) - x(:)'),[],2);
end













