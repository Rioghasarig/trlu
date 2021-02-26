function mcarrqr(m,n,blkMax,A,lda,jpvt,tau,info)
coder.ceval('ccarrqr', int32(m), int32(n), int32(blkMax), ...
             coder.ref(A), int32(lda), coder.ref(jpvt),...
             coder.ref(tau), int32(info));
         
end

