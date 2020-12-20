function [mylu] = factor(matfile)
load(matfile);
A = Problem.A;
[m,~] = size(A);
nrank = round(m/10); 
mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'nzinit', 50000000);
end

