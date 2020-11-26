%warning('off','all');
myDir = './test_data';
matFiles = dir(fullfile(myDir, '*.mat'));
nmat = length(matFiles);

% Table variables %
matname = strings(nmat,1); 
matsize = zeros(nmat,1);
facrank = zeros(nmat,1);
factime = zeros(nmat,1); 
facerror = -1*ones(nmat,1);
nswaps = zeros(nmat,1); 
swaperror = -1*ones(nmat,1);
swaptime = zeros(nmat,1); 
Ucond = zeros(nmat,1); 
swaps = 20; 
Unnz = zeros(nmat,1);
detA11 = zeros(nmat,swaps+1); 
for k = 1:1
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
    facrank(k) = nrank; 
    fprintf('Size: %d x %d\n', m, n);
    fprintf('Factorization Rank: %d\n', nrank); 
   %% Compute the matrix factorization
    % Finds permutations P, Q and matrices L and U such that PAQ  = L*U
    try
        tic; 
        mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'nzinit', 50000000);
        factime(k) = toc; 
    catch e
        fprintf('FACTORIATION FAILED: %s\n', e.message);
        continue
    end
    fprintf('Factorization Time: %.5f\n', factime(k));
    
    % Record the log of the determinant of the nrank x nrank upper-left
    % submatrix of PAQ = LU. To confirm the determinant is increasing
    detA11(k,1) = sum(log10(abs(mylu.diagU))); 
    
    enorm =mylu.facerror();
    facerror(k) = enorm;
    fprintf('Error: %.15f\n', enorm); 

    tic; 
    for i = 1:swaps
        
        [alpha, s_r, s_c] = mylu.maxS(); 
        [beta, a_r, a_c] = mylu.maxA11inv(alpha,s_r, s_c); 

        if alpha*beta < 2
           break
        end
        mylu.swapFac(a_r,a_c,s_r,s_c); 
        nswaps(k) = nswaps(k) + 1; 
        detA11(k,i+1)  = sum(log10(abs(mylu.diagU)));
    end
    swaptime(k) = toc; 
    fprintf('Number of Swaps: %d\n', nswaps(k)); 
    enorm = mylu.facerror(); 
    swaperror(k) = enorm;
    fprintf('Swap Error: %.15f\n', enorm); 
    fprintf('Swap TIme: %.5f\n', swaptime(k));
    
    % Plot swap figure
    plot(detA11(k,1:nswaps(k)+1));
    title(baseFileName, 'Interpreter', 'none');
    saveas(gcf, join(['figs/',baseFileName, '.fig']));
    
    diagU = mylu.diagU();
    Ucond(k) =max(abs(diagU))/min(abs(diagU));   
    fprintf('Ucond: %d\n', Ucond(k)); 
    Unnz(k) = nnz(mylu.U(1:nrank,:));
    fprintf('Unnz: %d\n', Unnz(k)); 
    
    T = table(matname, matsize, facrank,factime, facerror,nswaps,...
                    swaptime,swaperror, Ucond, Unnz,...
                    'VariableNames', {'Name', 'Size', 'Rank','FacTime', 'FacError', 'NumSwaps', ...
                   'SwapTIme','SwapError','Ucond','Unnz',});
writetable(T, 'results.csv'); 
end

