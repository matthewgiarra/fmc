MATLAB="/Applications/MATLAB_R2013a.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/matthewgiarra/.matlab/R2013a"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for testCompile" > testCompile_mex.mki
echo "CC=$CC" >> testCompile_mex.mki
echo "CFLAGS=$CFLAGS" >> testCompile_mex.mki
echo "CLIBS=$CLIBS" >> testCompile_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> testCompile_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> testCompile_mex.mki
echo "CXX=$CXX" >> testCompile_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> testCompile_mex.mki
echo "CXXLIBS=$CXXLIBS" >> testCompile_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> testCompile_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> testCompile_mex.mki
echo "LD=$LD" >> testCompile_mex.mki
echo "LDFLAGS=$LDFLAGS" >> testCompile_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> testCompile_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> testCompile_mex.mki
echo "Arch=$Arch" >> testCompile_mex.mki
echo OMPFLAGS= >> testCompile_mex.mki
echo OMPLINKFLAGS= >> testCompile_mex.mki
echo "EMC_COMPILER=" >> testCompile_mex.mki
echo "EMC_CONFIG=optim" >> testCompile_mex.mki
"/Applications/MATLAB_R2013a.app/bin/maci64/gmake" -B -f testCompile_mex.mk
