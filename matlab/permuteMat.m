function [A] = permuteMat(A,P1,Q1)
plen = length(P1);
qlen = length(Q1); 
A(1:plen,:) = A(P1, :);
A(:, 1:qlen) = A(:,Q1); 
end

