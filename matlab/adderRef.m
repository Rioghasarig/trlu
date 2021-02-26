function in1 = adderRef(m, n, blkMax, in1,lda, jpvt, in2, info)
% the input numel(in1) is converted to integer type 
% to match the cAdd function signature
coder.ceval('cAdd', m,n,blkMax,...
            coder.ref(in1), lda, coder.ref(jpvt),...
            coder.ref(in2), info );
end
