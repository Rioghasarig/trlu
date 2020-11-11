function [L U p q] = lucp(A,rank, tol, pm_opt)
%LUCP     LU factorization with complete pivoting.
%
% To compute the LU factorization under default settings:
%
%  [L U p q] = lucp(A)
%
% This produces a factorization such that L*U = A(p,q).  Vectors p and q
% permute the rows and columns, respectively.
%
% The pivot tolerance can be controlled:
%
%  [L U p q] = lucp(A,tol)
%
% The algorithm will terminate if the absolute value of the pivot is less
% than tol.
%
% Permutation matrices can be generated:
%
%  [L U P Q] = lucp(A,tol,'matrix')
%  [L U P Q] = lucp(A,tol,'sparse')
%
% The first will generate full permutation matrices P and Q such that 
% L*U = P*A*Q.  The second generates sparse P and Q.
%
% If A is sparse, L and U will be sparse.
%
% This function works on non-square matrices.
%
% Input:
%  A = matrix
%  tol = pivot tolerance (defualt is 1e-10)
%  pm_opt = permutation output options
%         = 'vector' for permutation vectors, L*U = A(p,q), defualt
%         = 'matrix' for full permutation matrices, L*U = P*A*Q
%         = 'sparse' for sparse permutation matrices, L*U = P*A*Q
%
% Output:
%  L = lower triangular factor
%  U = upper triangular factor
%  p = row permutation
%  q = column permutation
%
% Reference:
%  Algorithm 3.4.2, Matrix Computations, Golub and Van Loan. (3rd ed)
%  
% Other Implementations:
%  Gaussian Elimination using Complete Pivoting.  Alvaro Moraes.
%  http://www.mathworks.com/matlabcentral/fileexchange/25758
%
%  Gauss elimination with complete pivoting.  Nickolas Cheilakos.
%  http://www.mathworks.com/matlabcentral/fileexchange/13451
%  (Does not work with rectangular A)
%
%  Rank Revealing Code.  Leslie Foster.
%  http://www.math.sjsu.edu/~foster/rankrevealingcode.html
%  (Uses mex libraries for computation)
%


%
% 2010-03-28 (nwh) first version.
% 2010-04-14 (nwh) added more references.
% 2010-04-24 (nwh) perform final column swap so U is well conditioned
%

%
% License: see license.txt.
%

% handle optional inputs
if nargin < 3 || isempty(tol)
  tol = 1e-10;
end

if nargin < 4 || isempty(pm_opt)
  pm_opt = 'vector';
end

if strcmp(pm_opt,'vector')
  pm_flag = false;
  sp_flag = false;
elseif strcmp(pm_opt,'matrix')
  pm_flag = true;
  sp_flag = false;
elseif strcmp(pm_opt,'sparse')
  pm_flag = true;
  sp_flag = true;
else
  error('lucp:invalidOption','''%s'' is an un recognized option.',pm_opt)
end

[n m] = size(A);

% pivot vectors
p = (1:n)';
q = (1:m)';

% temp storage
rt = zeros(m,1); % row temp
ct = zeros(n,1); % col temp
t = 0; % scalar temp
maxnz = 10e7;
lv = zeros(maxnz,1);
li = zeros(maxnz,1);
lj = zeros(maxnz,1); 
lenl = 0; 

uv = zeros(maxnz,1);
ui = zeros(maxnz,1);
uj = zeros(maxnz,1); 
lenu = 0; 
S = A; 
%[ai, aj, av] = find(A); 
%A = sparse(ai,aj, av, n, m, 50000000);
for k = 1:rank
  % determine pivot
  %[cv ri] = max(abs(A(p(k:n),q(k:m))));
  srows = 1:n-(k-1);
  scols = 1:m-(k-1);
  [cv ri] = max(abs(S));
  [rv ci] = max(cv);
  %maxS = S(ri,ci); 
  rp = ri(ci)+k-1;
  cp = ci+k-1;
  rs = ri(ci); 
  cs = ci; 
  

  si = ri(ci); 
  % swap row
  t = p(k);
  p(k) = p(rp);
  p(rp) = t;
  %S([1,rs],:) = S([rs,1],:); 
  srows([rs,1]) = srows([1,rs]); 
  % swap col
  t = q(k);
  q(k) = q(cp);
  q(cp) = t;
  %S(:, [1,ci]) = S(:,[ci,1]); 
  scols([cs,1]) = scols([1,cs]); 

  
  if abs(S(rs,cs)) >= tol
   %% Add to u-file
    %[i, j, v] = find(A(p(k),q(k:m)));
    [i,j,v] = find(S(rs,scols)); 
    unnz = length(i);
    urng = lenu+1:lenu+unnz;
    ui(urng) = p(k);
    uj(urng) = q(j+k-1); 
    uv(urng) = v; 
    lenu = lenu + unnz; 
    
    %% Add to l-file
   % [i, j, v] = find(A(p(rows),q(k)));
    [i, j, v] = find(S(srows(2:end),cs));
    
    lnnz = length(i);
    lrng = lenl+1:lenl+lnnz;
    li(lrng) = p(i+k);
    lj(lrng) = p(k); 
    %lv(lrng) = v/A(p(k),q(k)); 
    lv(lrng) = v/S(rs,cs); 
    lenl = lenl+ lnnz; 
    Lelm = sparse(i,1,lv(lrng),n-k,1);     

    %% Update S
    %S = A(p(rows),q(cols))-Lelm*A(p(k),q(cols));
    S = S(srows(2:end),scols(2:end)) - Lelm*S(rs,scols(2:end)); 

  else
    % stop factorization if pivot is too small
    break
  end
end

% % final column swap if m > n
% if m > n
%   % determine col pivot
%   [cv ci] = max(abs(A(n,n:m)));
%   cp = ci+n-1;
%   
%   % swap col
%   t = q(n);
%   q(n) = q(cp);
%   q(cp) = t;
%   ct = A(:,n);
%   A(:,n) = A(:,cp);
%   A(:,cp) = ct;
% end

% produce L and U matrices
% these are sparse if L and U are sparse
%l = min(n,m);
%l = rank;

L =sparse(li(1:lenl),lj(1:lenl),lv(1:lenl), n, n) + speye(n,n);
%U = triu(A(p(1:n),q(1:m)));
U = sparse(ui(1:lenu),uj(1:lenu),uv(1:lenu),m,n); 
% produce sparse permutation matrices if desired
if pm_flag
  p = sparse(1:n,p,1);
  q = sparse(q,1:m,1);
end

% produce full permutation matrices if desired
if ~sp_flag
  p = full(p);
  q = full(q);
end
