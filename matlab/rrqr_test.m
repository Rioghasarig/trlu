
loadlibrary(' /global/home/groups/consultsw/sl-7.x86_64/modules/lapack/3.8.0/libblas.so');
A = gallery('dramadah', 10);
m = int32(10);
n = int32(10);
blkMax = int32(10);
lda = int32(10);
jpvt = int32(zeros(10,1));
tau = zeros(10,1);
info = int32(0);

a_ptr = libpointer('doublePtr', A);
m_ptr = libpointer('int32Ptr', m);
n_ptr = libpointer('int32Ptr',n);
blkMax_ptr = libpointer('int32Ptr',blkMax);
lda_ptr = libpointer('int32Ptr',lda);
jpvt_ptr = libpointer('int32Ptr',jpvt);
tau_ptr = libpointer('doublePtr', tau);
info_ptr = libpointer('int32Ptr',info);

%mcarrqr_mex(m,n,blkMax,A,lda,jpvt,tau,info);
loadlibrary('libclusol',@libclusol_proto_glnxa64);
calllib('libclusol','ccarrqr',...
    m_ptr,...
    n_ptr,...
    blkMax_ptr,...
    a_ptr,...
    lda_ptr,...
    jpvt_ptr,...
    tau_ptr,...
    info_ptr);