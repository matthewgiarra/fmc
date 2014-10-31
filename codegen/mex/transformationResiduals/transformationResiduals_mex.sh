MATLAB="/usr/local/MATLAB/R2013a"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/home/voodoo/.matlab/R2013a"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for transformationResiduals" > transformationResiduals_mex.mki
echo "CC=$CC" >> transformationResiduals_mex.mki
echo "CFLAGS=$CFLAGS" >> transformationResiduals_mex.mki
echo "CLIBS=$CLIBS" >> transformationResiduals_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> transformationResiduals_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> transformationResiduals_mex.mki
echo "CXX=$CXX" >> transformationResiduals_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> transformationResiduals_mex.mki
echo "CXXLIBS=$CXXLIBS" >> transformationResiduals_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> transformationResiduals_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> transformationResiduals_mex.mki
echo "LD=$LD" >> transformationResiduals_mex.mki
echo "LDFLAGS=$LDFLAGS" >> transformationResiduals_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> transformationResiduals_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> transformationResiduals_mex.mki
echo "Arch=$Arch" >> transformationResiduals_mex.mki
echo OMPFLAGS= >> transformationResiduals_mex.mki
echo OMPLINKFLAGS= >> transformationResiduals_mex.mki
echo "EMC_COMPILER=" >> transformationResiduals_mex.mki
echo "EMC_CONFIG=optim" >> transformationResiduals_mex.mki
"/usr/local/MATLAB/R2013a/bin/glnxa64/gmake" -B -f transformationResiduals_mex.mk
