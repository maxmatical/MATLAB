function [B] = Fast_MDS_B(E, F)
[m n_m] = size(F);
E_2 = E.^2;
B = zeros(m, n_m);
for i = 1:m
    for j = 1:n_m
        Ej_sq = sum(E_2, 2);
        B(i,j) = -1/2*(F(i,j)^2 - 1/m*Ej_sq(i));
    end
end

end