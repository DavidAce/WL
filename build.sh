#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-c] [-h ] [-j <num_threads>] [-l] [-m <mode>] [-t <target>]

-c               : Clear CMake files before build (delete ./build)
-h               : Help. Shows this text.
-j <num_threads> : Number of threads used by CMake (default $(nproc))
-l               : Clear downloaded libraries before build (i.e. delete ./libs)
-i               : Set install prefix (default: install)
-m <mode>        : Release   | Debug | (default = Release)
-t <target>      : DMRG++    | all   | any test target | (default = all)
EOF
  exit 1
}


target="all"
mode="Release"
clear_cmake=""
clear_libs=""
threads=$(nproc)
install="install"

while getopts chj:li:m:t: o; do
    case $o in
        (c) clear_cmake="true";;
        (h) usage ;;
        (j) threads=$OPTARG;;
        (l) clear_libs="true";;
        (i) install=$OPTARG;;
        (m) mode=$OPTARG;;
        (t) target=$OPTARG;;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage ;;
  esac
done
shift "$((OPTIND - 1))"


if [ "$clear_cmake" = "true" ]
then
    echo "Clearing CMake files from build."
	rm -f ./build/$mode/CMakeCache.txt
fi

echo "Starting Build"
echo "Target          :   $target"
echo "Build threads   :   $threads"
echo "Mode            :   $mode"
echo "Install prefix  :   $install"


cmake -E make_directory build/$mode
cd build/$mode
cmake -DCMAKE_BUILD_TYPE=$mode -DCMAKE_INSTALL_PREFIX=$install ../../
cmake --build . --target $target --parallel $threads
