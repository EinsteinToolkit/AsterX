#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors

################################################################################
# Search
################################################################################

if [ -z "${KADATH_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "KADATH selected, but KADATH_DIR not set."
    
    # Check for a previously built FUKA library
    if [ -n "${HOME_KADATH}" ] && [ -e "${HOME_KADATH}/lib/libkadath.a" ]; then
        echo "Installation of Kadath found at ${HOME_KADATH}"
        echo "Setting KADATH_DIR = $HOME_KADATH"
        echo "To force compiling the FUKA library, set KADATH_DIR = BUILD"
        KADATH_DIR=${HOME_KADATH}
    fi
    echo "END MESSAGE"
else
    echo "BEGIN MESSAGE"
    echo "Using KADATH in ${KADATH_DIR}"
    echo "END MESSAGE"
fi

THORN=KadathThorn

################################################################################
# Build
################################################################################

if [ -z "${KADATH_DIR}"                                                 \
     -o "$(echo "${KADATH_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]
then
    echo "BEGIN MESSAGE"
    echo "Building Frankfurt University/KADATH from git repo..."

    # Set locations
    SRCDIR="$(dirname $0)"
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${KADATH_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        INSTALL_DIR=${KADATH_INSTALL_DIR}
    fi
    echo "Installing Frankfurt University/KADATH into ${INSTALL_DIR}"
    echo "END MESSAGE"
    KADATH_BUILD=1
    KADATH_DIR=${INSTALL_DIR}
else
    KADATH_BUILD=
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    if [ ! -e ${DONE_FILE} ]; then
        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
        date > ${DONE_FILE}
    fi
fi

################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "KADATH_BUILD          = ${KADATH_BUILD}"
echo "KADATH_INSTALL_DIR    = ${KADATH_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set options
KADATH_INC_DIRS="${KADATH_DIR}/include ${KADATH_DIR}/include/Kadath_point_h"
KADATH_LIB_DIRS="${KADATH_DIR}/lib"
KADATH_LIBS="kadath"

echo "BEGIN MESSAGE"
echo "KadathThorn: Detect settings - KADATH_DIR : ${KADATH_DIR}"
echo "KadathThorn: Detect settings - KADATH_LIB_DIRS : ${KADATH_LIB_DIRS}"
echo "KadathThorn: Detect settings - KADATH_LIBS : ${KADATH_LIBS}"
echo "END MESSAGE"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "KADATH_DIR      = ${KADATH_DIR}"
echo "KADATH_INC_DIRS = ${KADATH_INC_DIRS}"
echo "KADATH_LIB_DIRS = ${KADATH_LIB_DIRS}"
echo "KADATH_LIBS     = ${KADATH_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(KADATH_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(KADATH_LIB_DIRS)'
echo 'LIBRARY           $(KADATH_LIBS)'
