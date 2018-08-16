function [x y] = Nys(E, F, d)
A = Fast_MDS_A(E);
[u gamma v] = svd(A);
ud = u(:, 1:d);
B = Fast_MDS_B(E, F);
gamma = gamma(1:d, 1:d);
x =  sqrt(gamma)*ud';
y = pinv(x')*B;
end
