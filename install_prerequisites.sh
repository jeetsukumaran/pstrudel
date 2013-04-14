#!/bin/sh
# set -x

if test "$1" == "--force"
then
    echo "Forcing rebuilding of all dependencies"
    unset NCL_ROOT
    unset NCL_PREFIX
fi

env_filename="$PWD/pstrudel_deps_env.sh"
echo '#!/bin/sh' > "${env_filename}"
PSTRUDEL_ROOT="$PWD"
echo "export PSTRUDEL_ROOT=${PSTRUDEL_ROOT}" >> "${env_filename}"

LIB_ROOT="$PSTRUDEL_ROOT/ext"
echo "export LIB_ROOT=${LIB_ROOT}" >> "${env_filename}"
if ! test -d $LIB_ROOT
then
    mkdir $LIB_ROOT
fi

echo $OSTYPE | grep darwin >/dev/null
is_linux=$?

export MAKE="make"

################################################################################
# NCL
################################################################################
if test -z $NCL_PREFIX
then
    if test -z $NCL_ROOT
    then
        cd $LIB_ROOT
        if ! test -d v2.1
        then
            svn checkout https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1 || exit 1
        fi
        NCL_ROOT="$LIB_ROOT/v2.1"
    fi
    echo "export NCL_ROOT=${NCL_ROOT}" >> "${env_filename}"
    cd "$NCL_ROOT" || exit 1
    mkdir build
    sh "./bootstrap.sh" || exit 1
    cd build
    export NCL_PREFIX="${NCL_ROOT}/build/installed"
    CXXFLAGS="-std=c++11" ../configure --prefix="${NCL_PREFIX}" || exit 1
    CXXFLAGS="-std=c++11" ${MAKE} || exit 1
    ${MAKE} install || exit 1
    ${MAKE} install check || exit 1
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    fi
    ${MAKE} check || exit
    cd $PSTRUDEL_ROOT
else
    if test ${is_linux} -eq 0
    then
        export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"
    fi
fi
echo "export NCL_PREFIX=${NCL_PREFIX}" >> "${env_filename}"
if test ${is_linux} -eq 0
then
    echo 'export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"'>> "${env_filename}"
else
    echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${NCL_PREFIX}/lib/ncl"'>> "${env_filename}"
fi

################################################################################
# Convenience
################################################################################

cat << 'EOF' > cfgcommand.sh
#! /bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if test -z $NCL_DIR
then
    source ${SCRIPT_DIR}/pstrudel_deps_env.sh
fi
if [[ $1 == "debug" ]]
then
    CXXFLAGS="-g -O0 -Wall -DPSTRUDEL_DEBUG_PRINTING -DDEBUGGING_MODE" ${SCRIPT_DIR}/configure --prefix=`pwd`/installed "--with-ncl=${NCL_PREFIX}"
    echo
    echo ---
    echo Debugging configuration complete.
elif [[ $1 == "profile" ]]
then
    CXXFLAGS="-O3 -pg -g -Wall" ${SCRIPT_DIR}/configure --prefix=`pwd`/installed "--with-ncl=${NCL_PREFIX}"
    echo
    echo ---
    echo Profile configuration complete.
else
    CXXFLAGS="-O3 -Wall" ${SCRIPT_DIR}/configure --prefix=`pwd`/installed "--with-ncl=${NCL_PREFIX}"
    echo
    echo ---
    echo Release configuration complete.
fi
EOF
chmod a+x cfgcommand.sh

echo
echo ---
echo "Run './bootstrap.sh' first if this is a clean checkout (or if you do not yet have a './configure' file in the project root directory."
echo "Configure by invoking 'cfgcommand.sh debug', 'cfgcommand.sh profile', or 'cfgcommand.sh release' (script should be in the project root directory, but invoked from the build directory)."

