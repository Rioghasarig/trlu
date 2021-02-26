MATLAB="/global/software/sl-7.x86_64/modules/tools/matlab/r2017b"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/global/home/users/oekenta/.matlab/R2017b"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for mcarrqr" > mcarrqr_mex.mki
echo "CC=$CC" >> mcarrqr_mex.mki
echo "CFLAGS=$CFLAGS" >> mcarrqr_mex.mki
echo "CLIBS=$CLIBS" >> mcarrqr_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> mcarrqr_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> mcarrqr_mex.mki
echo "CXX=$CXX" >> mcarrqr_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> mcarrqr_mex.mki
echo "CXXLIBS=$CXXLIBS" >> mcarrqr_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> mcarrqr_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> mcarrqr_mex.mki
echo "LDFLAGS=$LDFLAGS" >> mcarrqr_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> mcarrqr_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> mcarrqr_mex.mki
echo "Arch=$Arch" >> mcarrqr_mex.mki
echo "LD=$LD" >> mcarrqr_mex.mki
echo OMPFLAGS= >> mcarrqr_mex.mki
echo OMPLINKFLAGS= >> mcarrqr_mex.mki
echo "EMC_COMPILER=gcc" >> mcarrqr_mex.mki
echo "EMC_CONFIG=debug" >> mcarrqr_mex.mki
"/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/glnxa64/gmake" -B -f mcarrqr_mex.mk
