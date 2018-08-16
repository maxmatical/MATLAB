function [W H] = naiveNMF(A)
% initializing values
[m,n] = size(A);
k = min(n,m);
sigma = zeros(m);
U = zeros(m,m);
V = zeros(n,m);
% NMF based on power method
for mu = 1:k
    % generating variables for power method
    u = normrnd(0, 1, m, 1); 
    tolerance = 100;
    s = 0;
    v = zeros(m,1);
    % power method
    while tolerance >= 0.00001
        u_old = u;
        s_old = s;
        v_old = v;
        v = A'*u/norm(A'*u);
        s = norm(A*v);
        u = A*v/norm(A*v);
        tolerance = max(abs(s-s_old)); %,abs(u-u_old) , abs(v-v_old));
    end
    A = A - u*s*v';
    % replacing negative entries in A with 0
    for i = 1:m
        for j = 1:n
            if A(i,j) <0
                A(i,j)=0;
            end
        end
    end
    U(:,mu) = u;
    V(:,mu) = v;
    sigma(mu,mu) = s;    
end
W = U;
H = sigma*V';
end