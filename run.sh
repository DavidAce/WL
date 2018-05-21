#!/bin/bash
PROGNAME=$0

usage() {
  cat << EOF >&2

Usage            : $PROGNAME [-a <arg>] [-h] [-m <mode>] [-n] [-t <target>] [-v]

-a <arg>         : Input argument to WL (not implemented)
-h               : Help. Shows this text.
-m <mode>        : Release   | Debug  | (default = Release) -- (Use the same mode you used with build.sh)
-n               : Set number of threads to run with MPI
-t <target>      : WL    | all | (default = WL)
-v               : Use valgrind to find memory leaks
EOF
  exit 1
}

target="WL"
mode="Release"
arg=""
numcores="4"
valgrind=""
while getopts a:hm:n:t:v o; do
      case $o in
        (a) arg=$OPTARG;;
        (h) usage ;;
        (m) mode=$OPTARG;;
        (n) numcores=$OPTARG;;
        (t) target=$OPTARG;;
        (v) valgrind="valgrind --tool=memcheck --leak-check=full -v";;
        (:) echo "Option -$OPTARG requires an argument." >&2 ; exit 1 ;;
        (*) usage
  esac
done
shift "$((OPTIND - 1))"

bindtocore="-bind-to core:overload-allowed"
if [[ "$OSTYPE" == "darwin"* ]]; then
    bindtocore=""
fi

echo "Running command:  $valgrind mpirun -n $numcores $bindtocore ./build/$mode/$target $arg"
ulimit -c unlimited
$valgrind mpirun -n $numcores -bind-to core:overload-allowed ./build/$mode/$target $arg
