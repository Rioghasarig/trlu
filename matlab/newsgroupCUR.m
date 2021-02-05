load('20newsgroup-py.mat');
vocab_data = fileread('vocab-py.txt');
vocab = strsplit(vocab_data); 

A = Problem.train_data;
%A = (A~=0);
train_labels = Problem.train_labels;
c1 = find(train_labels == 0);
c2 = find(train_labels == 1); 
A1 = A(c1,:);
A2 = A(c2,:);

swaps = 30;
lmaxrs = zeros(swaps,1);
lmaxcs = zeros(swaps,1);
f = 1.1; 
nrank = 30;
mylu1 = lusol_obj(A1, 'pivot', 'TCP', 'rank', nrank, 'Ltol1', 10, 'Ltol2', 10, 'nzinit', 50000000);
mylu2 = lusol_obj(A2, 'pivot', 'TCP', 'rank', nrank, 'Ltol1', 10, 'Ltol2', 10, 'nzinit', 50000000);


prewordlist1 = {vocab{mylu1.aq(1:nrank)}};
prewordlist2 = {vocab{mylu2.aq(1:nrank)}};
rng(5);    
nswaps = 0;
for i = 1:swaps
    fprintf('Swap: %d\n', i); 
    [alpha, s_r, s_c] = mylu1.maxS(); 
    [beta, a_r, a_c] = mylu1.maxA11inv(alpha,s_r, s_c);
    
    fprintf('alpha*beta: %.15f\n', abs(alpha*beta));
    if abs(alpha*beta) < f
       break
    end
    nswaps = nswaps+1;
    [lmaxr,lmaxc] = mylu1.swapFac(a_r,a_c,s_r,s_c);
    lmaxrs(nswaps) = lmaxr;
    lmaxcs(nswaps) = lmaxc;
end


rng(5);    
for i = 1:swaps
    fprintf('Swap: %d\n', i); 
    [alpha, s_r, s_c] = mylu2.maxS(); 
    [beta, a_r, a_c] = mylu2.maxA11inv(alpha,s_r, s_c);

    L = mylu2.mulL(eye(584));
    U = mylu2.getU();
    E = mylu.A - L*U;
    S = mylu.A(nrank+1:end,nrank+1:end) - mylu.mulA22(eye(n-nrank));
    alphaT = max(max(abs(S)));
    A11 = full(mylu.A([1:nrank, nrank+s_r], [1:nrank, nrank+s_c]));
    betaT = max(max(abs(A11^-1)));   
    
    
    fprintf('alpha, alphaT: %.5f , %.5f\n', alpha, alphaT);
    fprintf('beta, betaT: %.5 , %.5f\n', beta, betaT); 
    if abs(alpha*beta) < f
       break
    end      

    mylu2.swapFac(a_r,a_c,s_r,s_c);
end

wordlist1 = {vocab{mylu1.aq(1:nrank)}};
wordlist2 = {vocab{mylu2.aq(1:nrank)}};

B = Problem.test_data;
test_labels = Problem.test_labels;

B1 = B(:,mylu1.aq);
B1norm = vecnorm(B1,2,2);
X1 = B1' ./B1norm'; 
fprintf('Computing B1 Error\n');
[B1x,inform,resid1] = mylu1.solveUt(X1);


B2 = B(:,mylu2.aq);
B2norm = vecnorm(B2,2,2);
X2 = B2' ./B2norm';
fprintf('Computing B2 Error\n');
[B2x,inform,resid2] = mylu2.solveUt(X2);


d1 = find(test_labels == 0);
d2 = find(test_labels == 1);

err1 = nnz(resid2(d1) >= resid1(d1))/length(d1);
err2 = nnz(resid1(d2) >= resid2(d2))/length(d2);
