function [A B] = Fast_MDS(D, m)
n = size(D, 2);
E = D(1:m, 1:m);
E_2 = E.^2;
F = D(1:m, (m+1):end);
A = zeros(m,m);
B = zeros(m, n-m);
for i = 1:m
    for j = 1:m
        Ej_sq = sum(E_2, 1);
        Ei_sq = sum(E_2, 2);
        Eij_sq = sum(sum(E_2,1),2);
        A(i,j) = -1/2*(E(i,j)^2- 1/m*Ei_sq(i) - 1/m*Ej_sq(j) + 1/m^2*Eij_sq);
        B(i,j) = -1/2*(F(i,j)^2 - 1/m*Ej_sq(j));
    end
end

end