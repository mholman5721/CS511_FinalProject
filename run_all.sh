#!/bin/bash

CUR_DIR=`pwd`
PY_OMP_SCRIPT=${CUR_DIR}/pythonHelperScripts/changeValsOmp.py
PY_SCRIPT=${CUR_DIR}/pythonHelperScripts/changeVals.py
RUN_DIRS=("gaussSBC_Test_ITERATIVE" "gaussSBC_Test_RING" "gaussSBC_TestOMP_RING" "jacobiIterative_Test" "jacobiBC_Test" "jacobiBC_TestOMP")
EXECUTABLES=("testGaussSBC" "testGaussSBC_Ring" "testGaussSBC_OMP" "testJacobiIterative" "testJacobiBC" "testJacobiBC_OMP")

DIMS_COUNT=10
THREADS_COUNT=1
MAX_THREADS=4
MAX_DIMS=300

CUR_EXECUTABLE=0

for dir in ${RUN_DIRS[@]};
do
    cd $dir
    echo $dir
    while [[ ${THREADS_COUNT} -le ${MAX_THREADS} ]]
    do
	  while [[ $DIMS_COUNT -le ${MAX_DIMS} ]]
	  do
	      if [[ $CUR_EXECUTABLE = 2 ]] || [[ $CUR_EXECUTABLE = 5 ]]; then
		 python ${PY_OMP_SCRIPT} ${CUR_DIR}/${dir}/const.h ${DIMS_COUNT} ${THREADS_COUNT}
	      else
		 python ${PY_SCRIPT} ${CUR_DIR}/${dir}/const.h ${DIMS_COUNT}
	      fi
	      
	      make
	      ./${EXECUTABLES[${CUR_EXECUTABLE}]}

	      ((DIMS_COUNT++))
	  done
	  mv "data_file_iterations.txt" $CUR_DIR/data_files/${EXECUTABLES[${CUR_EXECUTABLE}]}_${THREADS_COUNT}_iterations.txt
	  mv "data_file_time.txt" $CUR_DIR/data_files/${EXECUTABLES[${CUR_EXECUTABLE}]}_${THREADS_COUNT}_time.txt

	  ((THREADS_COUNT++))
	  DIMS_COUNT=1
    done
    ((CUR_EXECUTABLE++))
    THREADS_COUNT=1
    cd ${CUR_DIR}
done

    


