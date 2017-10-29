#!/usr/bin/env bash

OPTS=`getopt -o vhp:o:c: --long verbose,help,provider,mpich-dir:,ofi-dir:,csv-test-file: -n 'parse-options' -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

OFI_PROVIDER="sockets"
MPICH_DIR=""
OFI_DIR=""
PREV_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
CSV_TEST_FILE=""

while true; do
  case "$1" in
    -v | --verbose )  VERBOSE=true; shift ;;
    -h | --help )     HELP=true; shift ;;
    -p | --provider ) OFI_PROVIDER="$2"; shift 2 ;;
    --mpich-dir ) MPICH_DIR="$2"; shift 2 ;;
    -o | --ofi-dir ) OFI_DIR="$2"; shift 2 ;;
    -c | --csv-test-file ) CSV_TEST_FILE="$2"; shift 2 ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

export I_MPI_FABRICS=ofi
export FI_PROVIDER=$OFI_PROVIDER
echo FI_PROVIDER=$FI_PROVIDER

export LD_LIBRARY_PATH=$OFI_DIR/lib:$LD_LIBRARY_PATH

#echo "Building MPICH collective tests:"
#cd $MPICH_DIR/test/mpi/coll
#for i in *.c
#do
#  echo "Building $i >> ${i%.c}..."
#  mpiicc -I ../include/  -o "${i%.c}"  "$i" ../util/mtest.c ../util/mtest_datatype_gen.c ../util/mtest_datatype.c
#  echo "Running ${i%.c}..."
#  mpirun -n 4 "${i%.c}"
#done
#cd -

echo "MPICH tests:"
mpich_tests=()
cd $MPICH_DIR/test/mpi
while IFS="," read -r f1 f2 f3
do
  echo "Building $f1 >> ${f1%.c}..." && \
  mpiicc -I $MPICH_DIR/test/mpi/include/  -o "${f1%.c}" "$f1" $MPICH_DIR/test/mpi/util/mtest.c $MPICH_DIR/test/mpi/util/mtest_datatype_gen.c $MPICH_DIR/test/mpi/util/mtest_datatype.c && \
  mpich_tests=("${mpich_tests[@]}" "-n $f3 ${f1%.c} $f2")
done < "$CSV_TEST_FILE"

for i in "${mpich_tests[@]}";
do
  echo "Running mpirun $i";
  I_MPI_OFI_PROVIDER=$OFI_PROVIDER mpirun $i;
done

cd -

unset I_MPI_OFI_FABRICS
unset FI_PROVIDER
export LD_LIBRARY_PATH=$PREV_LD_LIBRARY_PATH

