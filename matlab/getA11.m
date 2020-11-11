function [A11] = getA11(trlu)

    p = trlu.p;
    q = trlu.q; 
    nrank = trlu.stats.nrank; 
    [m,n] = trlu.size(); 
    A11 = sparse(nrank, nrank); 
    for j = 1:nrank
        ej = zeros(n,1); 
        ej(q(j)) = 1;
        aj = trlu.mulA(ej); 
        A11(:,j) = aj(p(1:nrank));
    end

end
