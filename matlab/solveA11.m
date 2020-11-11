function [x] = solveA11(alu, b,nrank)
    [m,~] = size(alu); 
    p = alu.p;
    q = alu.q;
    U = alu.U;
 
  
    b2 = zeros(m, 1);
    b2(p(1:nrank)) = b;
    c = alu.solveL(b2);
    c2 = c(p(1:nrank));     
    U11 = U(p(1:nrank),q(1:nrank)); 

    x = U11\c2; 
end

