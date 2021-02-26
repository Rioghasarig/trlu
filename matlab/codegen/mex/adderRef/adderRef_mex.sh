MATLAB="/global/software/sl-7.x86_64/modules/tools/matlab/r2017b"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/global/home/users/oekenta/.matlab/R2017b"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for adderRef" > adderRef_mex.mki
echo "CC=$CC" >> adderRef_mex.mki
echo "CFLAGS=$CFLAGS" >> adderRef_mex.mki
echo "CLIBS=$CLIBS" >> adderRef_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> adderRef_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> adderRef_mex.mki
echo "CXX=$CXX" >> adderRef_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> adderRef_mex.mki
echo "CXXLIBS=$CXXLIBS" >> adderRef_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> adderRef_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> adderRef_mex.mki
echo "LDFLAGS=$LDFLAGS" >> adderRef_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> adderRef_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> adderRef_mex.mki
echo "Arch=$Arch" >> adderRef_mex.mki
echo "LD=$LD" >> adderRef_mex.mki
echo OMPFLAGS= >> adderRef_mex.mki
echo OMPLINKFLAGS= >> adderRef_mex.mki
echo "EMC_COMPILER=gcc" >> adderRef_mex.mki
echo "EMC_CONFIG=debug" >> adderRef_mex.mki
"/global/software/sl-7.x86_64/modules/tools/matlab/r2017b/bin/glnxa64/gmake" -B -f adderRef_mex.mk
