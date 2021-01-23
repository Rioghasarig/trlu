load('20newsgroup-py.mat');
vocab_data = fileread('vocab-py.txt');
vocab = strsplit(vocab_data); 

A = Problem.train_data;
%A = (A~=0);
train_labels = Problem.train_labels;
c1 = find(train_labels == 0);
c2 = find(train_labels == 1); 
train_labels = train_labels([c1,c2]);
A = A([c1,c2],:);
%Anorm = vecnorm(A,2,2);
%A = A ./ Anorm; 
swaps = 100;
lmaxrs = zeros(swaps,1);
lmaxcs = zeros(swaps,1);
f = 1.1; 
nrank = 1000;
mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'Ltol1', 10, 'Ltol2', 10, 'nzinit', 50000000);


prewordlist = {vocab{mylu.aq(1:nrank)}};
rng(5);    
nswaps = 0;
for i = 1:swaps
    fprintf('Swap: %d\n', i); 
    [alpha, s_r, s_c] = mylu.maxS2(); 
    [beta, a_r, a_c] = mylu.maxA11inv(alpha,s_r, s_c);
    
    fprintf('alpha*beta: %.15f\n', abs(alpha*beta));
    if abs(alpha*beta) < f
       break
    end
    nswaps = nswaps+1;
    [lmaxr,lmaxc] = mylu.swapFac(a_r,a_c,s_r,s_c);
    lmaxrs(nswaps) = lmaxr;
    lmaxcs(nswaps) = lmaxc;
    lenr = mylu.lenr_ptr.Value;
    lenu = mylu.stats.lenU;
    lrow = mylu.stats.lrow;
end


wordlist = {vocab{mylu.aq(1:nrank)}};
train_labels = train_labels(mylu.ap);
B = Problem.test_data;

test_labels = Problem.test_labels;

% B = B(:,mylu.aq);
% Bnorm = vecnorm(B,1,2);
% X = B' ./Bnorm'; 
% fprintf('Computing B1 Error\n');
% [Bx,inform,resid1] = mylu.solveUt(X);
% Bx = Bx(1:nrank,:); 
% ident = eye(2);
% e = ident(:,train_labels +1); 
% eBx = e*Bx; 
