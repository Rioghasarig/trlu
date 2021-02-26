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
detA = zeros(swaps,1);
lrs = zeros(swaps,3);
lcs = zeros(swaps,3);
errors= zeros(swaps,3);
f = 1.1; 
nrank = 1000;
mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'Ltol1', 10, 'Ltol2', 10, 'nzinit', 50000000);


prewordlist = {vocab{mylu.aq(1:nrank)}};
rng(5);    
nswaps = 0;
for i = 1:swaps
    fprintf('Swap: %d\n', i); 
    [alpha, s_r, s_c] = mylu.maxS(); 
    [beta, a_r, a_c] = mylu.maxA11inv(alpha,s_r, s_c);
    
%     diagU = mylu.diagU();
%     detA(i) = sum(log(abs(diagU))/log(f));   
%     L = mylu.L11;
%     U = mylu.getU();
% 
%     E0 = mylu.A(1:nrank,1:nrank) - L*U;
%     E1 = mylu.A(1:nrank,nrank+1:end) - L*mylu.U12;
%     E2 = mylu.A(nrank+1:end,1:nrank) - mylu.L21*U;
%     
%     e0 = max(max(abs(E0)));
%     e1 = max(max(abs(E1)));
%     e2 = max(max(abs(E2)));
%     errors(i,:) = [e0,e1,e2];
%     fprintf('e0: %.15f, e1: %.15f, e2: %.15f\n',...
%         full(e0), full(e1), full(e2)); 
%     fprintf('detA: %.15f\n', detA(i));
    if abs(alpha*beta) < f
       break
    end
    nswaps = nswaps+1;
    [lr, lc] = mylu.swapFac(a_r,a_c,s_r,s_c);
 %   fprintf('lr, lc: %d, %d\n', lr, lc);
    lrs(i,:) = lr;
    lcs(i,:) = lc;
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
