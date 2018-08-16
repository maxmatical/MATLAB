
load 'Simple.mat';

dist = zeros(6);
for i = 1:6
    for j = 1:6
        val = 0;
        for k = 1:3
            val = val + (X(k,i)-X(k,j))^2;
        end
        val = sqrt(val);
        dist (i,j) = val;
    end
end