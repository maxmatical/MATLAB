function [A] = Fast_MDS_A(E)
m = size(E,2);
A = zeros(m,m);
E_2 = E.^2;
for i = 1:m
    for j = 1:m
        Ej_sq = sum(E_2, 1);
        Ei_sq = sum(E_2, 2);
        Eij_sq = sum(sum(E_2,1),2);
        A(i,j) = -1/2*(E(i,j)^2- 1/m*Ei_sq(i) - 1/m*Ej_sq(j) + 1/m^2*Eij_sq);
    end
end
end