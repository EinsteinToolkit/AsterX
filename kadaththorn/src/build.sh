#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors

# SRCDIR = KadathThorn/src
SRCDIR="$(dirname $0)"
echo "KadathThorn - BUILD srcdir: ${SRCDIR}"
pushd ${SRCDIR}

BASEDIR_ABS=`cd ../; pwd -P`
pushd ${BASEDIR_ABS}
# Set locations
THORN=KadathThorn

NAME=fuka

# Setup temporary build and final installation directories
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${KADATH_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "Installing FUKA into ${KADATH_INSTALL_DIR}"
    INSTALL_DIR=${KADATH_INSTALL_DIR}
fi
echo "KadathThorn: Working in ${INSTALL_DIR}"

DONE_FILE=${SCRATCH_BUILD}/done/${THORN}

KADATH_DIR=${INSTALL_DIR}
if [ -e "${KADATH_DIR}/lib/libkadath.a" ]; then 
    echo "KadathThorn: FUKA has already been built."
    echo "KadathThorn: If you want to force rebuilding it, delete ${KADATH_DIR}/lib/libkadath.a"    
else

# Start from clean build and install directories
echo "KadathThorn: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "KadathThorn: Copying Frankfurt University/KADATH to ${BUILD_DIR} ..."
# Change to build directory while retaining the current directory under bash `dirs` list
pushd ${BUILD_DIR}

# Make a copy of the FUKA repo to the build location
cp -LR ${SRCDIR}/${NAME} ./

echo "KadathThorn: Configuring..."
# Change to directory containing the FUKA repo
cd ${NAME}

export HOME_KADATH=${BUILD_DIR}/${NAME}

# KADATH is built using cmake
# Local settings are stored in $HOME_KADATH/Cmake/CMakeLocal.cmake
# Here we generate our own CMakeLocal.cmake using build variables
# that should be populated by the user or the ETK builder
#
# If this is not robust, this can be commented out and the user set
# their own CMakeLocal.cmake file...or just build FUKA separately
rm ${HOME_KADATH}/Cmake/CMakeLocal.cmake
cat > ${HOME_KADATH}/Cmake/CMakeLocal.cmake <<EOF
set (GSL_LIBRARIES ${GSL_LIBS})
set (SCALAPACK_LIBRARIES ${LAPACK_LIBS})
set (FFTW_LIBRARIES ${FFTW3_LIBS})
set (BLAS_LIBRARIES ${BLAS_LIBS})

EOF

# Maybe this should only be for verbose/debug?
echo "KadathThorn: CMAKELIST is"
cat ${HOME_KADATH}/Cmake/CMakeLocal.cmake

# Move to build directory holding CMakeLists.txt
cd ${HOME_KADATH}/build_release
# Make build directory to store build files from CMake and make
mkdir -p build
cd build

# We only are concerned with importing initial data, not solving initial data
# Here we remove the do_newton from the list of compiled files to avoid
# all of its dependencies (SCALAPACK, etc)
sed -i -e  '/newton/d' ../CMakeLists.txt

# Generate Makefile based on CMakeLocal.cmake settings
cmake -DPAR_VERSION=ON -DCMAKE_BUILD_TYPE=Release ..

echo "KadathThorn: Building Frankfurt University/KADATH library..."
# Compile and build library $HOME_KADATH/lib/libkadath.a
${MAKE}

echo "KadathThorn: Installing Frankfurt University/KADATH library..."
mv ${BUILD_DIR}/${NAME}/lib ${INSTALL_DIR}
mv ${BUILD_DIR}/${NAME}/eos ${INSTALL_DIR}

echo "KadathThorn: Set environment variable HOME_KADATH to ${INSTALL_DIR} in job submissions"
echo "KadathThorn: To use FUKA's built-in equations of state"

# This could be useful later
# mv ${BUILD_DIR}/${NAME}/include        ${INSTALL_DIR}/

popd
echo "KadathThorn: Cleaning up..."
rm -rf ${BUILD_DIR}
fi

date > ${DONE_FILE}
echo "KadathThorn: Done."
