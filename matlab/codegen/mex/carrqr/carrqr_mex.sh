MATLAB="/global/software/sl-7.x86_64/modules/tools/matlab/r2017b"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/global/home/users/oekenta/.matlab/R2017b"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for carrqr" > carrqr_mex.mki
echo "CC=$CC" >> carrqr_mex.mki
echo "CFLAGS=$CFLAGS" >> carrqr_mex.mki
echo "CLIBS=$CLIBS" >> carrqr_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> carrqr_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> carrqr_mex.mki
echo "CXX=$CXX" >> carrqr_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> carrqr_mex.mki
echo "CXXLIBS=$CXXLIBS" >> carrqr_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> carrqr_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> carrqr_mex.mki
echo "LDFLAGS=$LDFLAGS" >> carrqr_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> carrqr_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> carrqr_mex.mki
echo "Arch=$Arch" >> carrqr_mex.mki
echo "LD=$LD" >> carrqr_mex.mki
echo OMPFLAGS= >> carrqr_mex.mki
echo OMPLINKFLAGS= >> carrqr_mex.mki
echo "EMC_COMPILER=gcc" >> carrqr_mex.mki
echo "EMC_CONFIG=optim" >> carrqr_mex.mki
"/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/glnxa64/gmake" -B -f carrqr_mex.mk
