CC="/usr/bin/gcc"
                CXX="g++"
                CFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -D_GNU_SOURCE -DMATLAB_MEX_FILE "
                CXXFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -D_GNU_SOURCE -DMATLAB_MEX_FILE "
                COPTIMFLAGS="-O -DNDEBUG"
                CXXOPTIMFLAGS="-O -DNDEBUG"
                CDEBUGFLAGS="-g"
                CXXDEBUGFLAGS="-g"
                LD="/usr/bin/gcc"
                LDXX="g++"
                LDFLAGS="-pthread -Wl,--no-undefined -Wl,-rpath-link,/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/glnxa64 -shared  -L"/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -Wl,--version-script,"/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/extern/lib/glnxa64/mexFunction.map""
                LDDEBUGFLAGS="-g"
