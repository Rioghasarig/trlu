%warning('off','all');
myDir = './small_mats';
matFiles = dir(fullfile(myDir, '*.mat'));
nmat = length(matFiles);

% Table variables %
matname = strings(nmat,1); 
matsize = zeros(nmat,1);
matnnz = zeros(nmat,1);
facrank = zeros(nmat,1);
factime = zeros(nmat,1); 
facerror = -1*ones(nmat,1);
facerror12 = -1*ones(nmat,1);
facerror21 = -1*ones(nmat,1); 
nswaps = zeros(nmat,1); 
swaperror = -1*ones(nmat,1);
swaptime = zeros(nmat,1); 
Ucond = zeros(nmat,2); 
swaps = 20; 
Unnz = zeros(nmat,1);
Lnnz = zeros(nmat, 1); 
detA11 = zeros(nmat,swaps+1); 
alphaError = zeros(nmat,swaps+1);
betaError = zeros(nmat,swaps+1); 
lmaxrs = zeros(nmat,swaps);
lmaxcs = zeros(nmat,swaps); 
f = 2;
for k = 6:6
    baseFileName = matFiles(k).name;
    nameLen = min(10, length(baseFileName));
    matname(k) = baseFileName(1:nameLen);
    fullFileName = fullfile(myDir,baseFileName);
    fprintf('\n=================\n')
    fprintf('Matrix: %s\n', baseFileName);
    load(fullFileName);
    A = Problem.A; 
    [m, n] = size(A);
    nrank = round(m/10);


    matsize(k) = m;
    matnnz(k) = nnz(A); 
    facrank(k) = nrank; 
    fprintf('Size: %d x %d\n', m, n);
    fprintf('Nnz: %d\n', matnnz(k)); 
    fprintf('Factorization Rank: %d\n', nrank); 
   %% Compute the matrix factorization
    % Finds permutations P, Q and matrices L and U such that PAQ  = L*U
    try
        tic; 
        mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'nzinit', 50000000);
        factime(k) = toc; 
        
            
        % Record the log of the determinant of the nrank x nrank upper-left
        % submatrix of PAQ = LU. To confirm the determinant is increasing 
        diagU = mylu.diagU();
        detA11(k,1) = sum(log10(abs(diagU))); 
        Ucond(k,1) =max(abs(diagU))/min(abs(diagU));   
        fprintf('Ucond1: %d\n', Ucond(k,1));

        enorm = mylu.facerror();
        enorm2 = mylu.facerror12();
        enorm3 = mylu.facerror21();
        facerror(k) = enorm;
        fprintf('Error: %.15f\n', enorm); 
        fprintf('Error2: %.15f\n', enorm2);
        fprintf('Error3: %.15f\n', enorm3); 
    catch e
        fprintf('FACTORIATION FAILED: %s\n', e.message);
        continue
    end
    fprintf('Factorization Time: %.5f\n', factime(k));


    try
        tic; 
        rng(5);    
        for i = 1:swaps
            fprintf('Swap: %d\n', i); 
            [alpha, s_r, s_c] = mylu.maxS(); 
            [beta, a_r, a_c] = mylu.maxA11inv(alpha,s_r, s_c); 
            
            
%             S = mylu.A(nrank+1:end,nrank+1:end) - mylu.L21*mylu.U12;
%             alphaT = full(max(max(abs(S))));          
%             A11 = mylu.A([1:nrank, nrank+s_r], [1:nrank, nrank+s_c]);
%             betaT = full(max(max(abs(A11^-1)))); 
%             
%             fprintf('alpha: %d , %d\n', alpha, alphaT);           
%             fprintf('beta: %d, %d\n', beta, betaT);
%             
%             alphaError(k,i) = abs(abs(alphaT)/abs(alpha));
%             betaError(k,i) = abs(abs(betaT)/abs(beta));             
            if abs(alpha*beta) < f
               break
            end            
            mylu.swapFac(a_r,a_c,s_r,s_c);

            detA11(k,i+1) = sum(log10(abs(mylu.diagU)));
            if detA11(k,i+1) - detA11(k,i) < log10(f)
                fprintf('WARNING: detA11 not increasing enough \n');
            end
            nswaps(k) = nswaps(k) + 1; 
        end
    catch e
        fprintf('Error during swaps: %s - %s\n', e.identifier, e.message);
        continue
    end
    
    try
       swaptime(k) = toc; 
       fprintf('Number of Swaps: %d\n', nswaps(k));
       %mylu.U = mylu.getU(); 
       %detA11(k,2) = sum(log10(abs(mylu.diagU))); 
       %fprintf('Det A11 diff: %.5f\n', detA11(k,2) - detA11(k,1)); 
       enorm = mylu.facerror(); 
       enorm2 = mylu.facerror12();
       enorm3 = mylu.facerror21();
       swaperror(k) = enorm;
       fprintf('Swap Error: %.15f\n', enorm); 
       fprintf('Swap Error2: %.15f\n', enorm2);
       fprintf('Swap Error3: %.15f\n', enorm3); 
       fprintf('Swap Time: %.5f\n', swaptime(k));
      
       diagU = mylu.diagU();
       Ucond(k,2) =max(abs(diagU))/min(abs(diagU));   
       fprintf('Ucond2: %d\n', Ucond(k,2)); 
       Unnz(k) = nnz(mylu.getU()); 
       fprintf('Unnz: %d\n', Unnz(k)); 
       Lnnz(k) = mylu.stats.lenL; 
       fprintf('Lnnz: %d\n', Lnnz(k)); 
       

       
       T = table(matname, matsize, matnnz, facrank,factime, facerror,nswaps,...
                       swaptime,swaperror, Ucond(:,1),Ucond(:,2), Unnz,Lnnz,...
                       'VariableNames', {'Name', 'Size', 'Matnnz', 'Rank','FacTime', 'FacError', 'NumSwaps', ...
                      'SwapTime','SwapError','Ucond1','Ucond2','Unnz', 'Lnnz'});
       writetable(T, 'results.csv'); 
    catch e
        fprintf('Error in post-swap process: %s', e.message);
        continue
    end
    
end

% alphaHist = [];
% j = 1;
% 
% for k = 1:nmat
%     for i = 1:nswaps(k)
%         alphaHist(j) = alphaError(k,i);
%         j = j+1;
%     end
% end
% 
% betaHist = [];
% j = 1; 
% for k = 1:nmat
%     for i = 1:nswaps(k)+1
%         betaHist(j) = betaError(k,i);
%         j = j + 1;
%     end
% end

