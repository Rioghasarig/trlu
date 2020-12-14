load('test_data/CoupCons3D.mat');
A = Problem.A;
[m,n] = size(A);
nrank = round(m/10);
mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'nzinit', 50000000);
mylu.swapRows(9354,nrank+42924);
