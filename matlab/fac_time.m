%warning('off','all');
myDir = './test_data';
matFiles = dir(fullfile(myDir, '*.mat'));
nmat = length(matFiles);

% Table variables %
matname = strings(nmat,1); 
matsize = zeros(nmat,1);
matnnz = zeros(nmat,1);
factime = zeros(nmat,1); 
facerror = -1*ones(nmat,1);
nswaps = zeros(nmat,1); 
swaperror = -1*ones(nmat,1);
swaptime = zeros(nmat,1); 
Ucond = zeros(nmat,2); 
swaps = 20; 
Unnz = zeros(nmat,1);
Lnnz = zeros(nmat, 1); 
detA11 = zeros(nmat,swaps+1); 
for k = 1:nmat
    baseFileName = matFiles(k).name;
    nameLen = min(10, length(baseFileName));
    matname(k) = baseFileName(1:nameLen);
    fullFileName = fullfile(myDir,baseFileName);
    load(fullFileName);
    A = Problem.A; 
    
    if size(A, 1) > 2*10^3
        continue
    end
    fprintf('\n=================\n')
    fprintf('Matrix: %s\n', baseFileName);

    [m, n] = size(A);
    

    matsize(k) = m;
    matnnz(k) = nnz(A); 
    fprintf('Size: %d x %d\n', m, n);
    fprintf('Nnz: %d\n', matnnz(k)); 
   %% Compute the matrix factorization
    % Finds permutations P, Q and matrices L and U such that PAQ  = L*U
    try
        tic; 
        mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', m, 'nzinit', 50000000);
        factime(k) = toc; 
        
            
        mylu.U = mylu.getU(); 
        diagU = mylu.diagU();
        detA11(k,1) = sum(log10(abs(diagU))); 
        Ucond(k,1) =max(abs(diagU))/min(abs(diagU));   
        fprintf('Ucond1: %d\n', Ucond(k,1));

        enorm =mylu.facerror();
        facerror(k) = enorm;
        fprintf('Error: %.15f\n', enorm); 
    catch e
        fprintf('FACTORIATION FAILED: %s\n', e.message);
        continue
    end
    fprintf('LUSOL Time: %.5f\n', factime(k));

    try
        tic
        [L, U, P,Q] = lu(A);
        ml_time = toc;
        fprintf('MATLAB TIME: %.5f\n', ml_time)
    end
        
    try

 
        T = table(matname, matsize, matnnz,factime,...
                        'VariableNames', {'Name', 'Size', 'Matnnz', 'FacTime'});
        writetable(T, 'fac_time.csv'); 
    catch e
        fprintf('Error in post-swap process: %s', e.message);
        continue
    end
end

