%warning('off','all');
myDir = './snap_mats';
matFiles = dir(fullfile(myDir, '*.mat'));
nmat = length(matFiles);

% Table variables %
matname = strings(nmat,1); 
matsize = zeros(nmat,1);
matnnz = zeros(nmat,1);
nswaps = zeros(nmat, 1); 

facrank = zeros(nmat,10);
factime = zeros(nmat,10);
L2error = -1*ones(nmat,10);
swaps = 20; 
f = 10;

for k = 1:nmat
    baseFileName = matFiles(k).name;
    nameLen = min(10, length(baseFileName));
    matname(k) = baseFileName(1:nameLen);
    fullFileName = fullfile(myDir,baseFileName);
    fprintf('\n=================\n')
    fprintf('Matrix: %s\n', baseFileName);
    
    load(fullFileName);% loads 'A' matrix
    %A = Problem.A;
    [m, n] = size(A);
    matsize(k) = m;
    matnnz(k) = nnz(A); 
    fprintf('Size: %d x %d\n', m, n);
    fprintf('Nnz: %d\n', matnnz(k)); 
    j = 0;
    for nrank = 100:100:1000
        j = j+1;
        %% Compute the factorization
        fprintf('Rank: %d\n', nrank);
        tic
        try
            mylu = lusol_obj(A, 'pivot', 'TCP', 'rank', nrank, 'nzinit', 50000000); 
        catch e
            fprintf('FACTORIATION FAILED: %s\n', e.message);
            break
        end

            
        try
            rng(5);    
            for i = 1:swaps
                fprintf('Rank: %d, Swap: %d\n', nrank, i); 
                [alpha, s_r, s_c] = mylu.maxS2(); 
                [beta, a_r, a_c] = mylu.maxA11inv(alpha,s_r, s_c); 
                fprintf('alpha*beta: %.15f\n', abs(alpha*beta));
                if abs(alpha*beta) < f
                   break
                end            
                mylu.swapFac(a_r,a_c,s_r,s_c);
                nswaps(k) = nswaps(k) + 1;
            end
        catch e
            fprintf('Error during swaps: %s - %s\n', e.identifier, e.message);
            continue
        end
        enorm = mylu.facerror();
        enorm2 = mylu.facerror12();
        enorm3 = mylu.facerror21();
        fprintf('Error: %.15f\n', enorm); 
        fprintf('Error2: %.15f\n', enorm2);
        fprintf('Error3: %.15f\n', enorm3); 
        
        factime(k,j) = toc; 
        fprintf('Rank: %d - Computing Error\n', nrank);
        L2error(k,j) = mylu.L2error();
        fprintf('Number of Swaps: %d\n', nswaps(k));
    end
end

