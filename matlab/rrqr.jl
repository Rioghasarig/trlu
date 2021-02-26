
using Libdl
using LinearAlgebra
using MAT
hlapack = dlopen("/global/home/groups/consultsw/sl-7.x86_64/modules/lapack/3.8.0/liblapack.so",RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
hblas = dlopen("/global/home/groups/consultsw/sl-7.x86_64/modules/lapack/3.8.0/libblas.so",RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
hrrqr = dlopen("/global/home/users/oekenta/trlu/matlab/rrqr.so")
homp = dlopen("/global/software/sl-7.x86_64/modules/langs/intel/2016.4.072/lib/intel64_lin/libiomp5",RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
trqrcp = dlsym(hrrqr,"trqrcp_")
io = open("facerror.out", "w")

qr_error = dlsym(hrrqr, "qr_error_")
mats = readdir("small_mats")
nmats = length(mats)
errTable = Dict()
for i=1:nmats
    file = matopen(string("small_mats/", mats[i]))
    Problem = read(file, "Problem")
    close(file)
  
    A = Matrix(Problem["A"])
    
    m = Ref{Int32}(size(A,1))
    n = Ref{Int32}(size(A,2))
    ldA = Ref{Int32}(size(A,1))
    jpvt = zeros(Int32,size(A,2),1)
    Y = zeros(size(A,1),size(A,2))
    ldY = Ref{Int32}(size(A,1))
    tau = zeros(size(A,2),1)
    R = zeros(size(A,1),size(A,2))
    ldR = Ref{Int32}(size(A,1))
    iseed = Int32[1,1,1,1]
    info = Ref{Int32}(0)

    facerrors = zeros(5,1)

    print("\n")
    println(i, ": ", mats[i])
    println("Size: ", size(A))

    write(io,"\n")
    write(io, string(i, ": ", mats[i], "\n"))
    write(io, string(mats[i],"\n"))
    write(io, string("Size: ", size(A), "\n"))
    for j=1:5
        facrank = 200*j
        k = Ref{Int32}(facrank)
        blkMax = Ref{Int32}(facrank)
        ccall(trqrcp,Cvoid, 
              (Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},Ptr{Float64},
               Ref{Int32},Ptr{Int32},Ptr{Float64},Ref{Int32}, Ptr{Float64},
               Ptr{Float64}, Ref{Int32},Ptr{Int32},Ref{Int32}),
              m,n,k,blkMax,A,ldA,jpvt,Y,ldY,tau,R,ldR,iseed,info)


        l = Ref{Int32}(size(A,1))
        err = ccall(qr_error,Float64,
                    (Ref{Int32},Ref{Int32},Ref{Int32},Ref{Int32},
                     Ptr{Float64}, Ref{Int32}, Ptr{Float64}, Ref{Int32},
                     Ptr{Float64}, Ref{Int32}),
                    m,n,k,l,A,ldA,Y,ldY,R,ldR)
        err = err*norm(R)/norm(A)
        println(facrank,": ",  err)
        write(io, string(facrank, ": ", err, "\n"))
        facerrors[j] = err

    end

    errTable[mats[i]] = facerrors
end
close(io)
