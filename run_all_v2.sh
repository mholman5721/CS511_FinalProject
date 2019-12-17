#!/bin/bash

CUR_DIR=`pwd`
PY_OMP_SCRIPT=${CUR_DIR}/pythonHelperScripts/changeValsOmp.py
PY_SCRIPT=${CUR_DIR}/pythonHelperScripts/changeVals.py
RUN_DIRS=("gaussSBC_TestOMP_ITERATIVE")
EXECUTABLES=("testGaussSBC_OMP_ITERATIVE")
#RUN_DIRS=("gaussSBC_Test_RING" "gaussSBC_Test_ITERATIVE") "gaussSBC_TestOMP_RING" 
#EXECUTABLES=("testGaussSBC_Ring" "testGaussSBC_Iterative") "testGaussSBC_OMP_RING" 

DIMS_COUNT=10
THREADS_COUNT=2
MAX_THREADS=2
MAX_DIMS=100

CUR_EXECUTABLE=0

for dir in ${RUN_DIRS[@]};
do
    cd $dir
    echo $dir
    while [[ ${THREADS_COUNT} -le ${MAX_THREADS} ]]
    do
	  while [[ $DIMS_COUNT -le ${MAX_DIMS} ]]
	  do
		python ${PY_OMP_SCRIPT} ${CUR_DIR}/${dir}/const.h ${DIMS_COUNT} ${THREADS_COUNT}
		#python ${PY_SCRIPT} ${CUR_DIR}/${dir}/const.h ${DIMS_COUNT} ${THREADS_COUNT}

		make clean
		sleep 0.1
	    make
		sleep 0.1
	    ./${EXECUTABLES[${CUR_EXECUTABLE}]}
		sleep 0.1
		
	    ((DIMS_COUNT++))
	  done
	  mv "data_file_iterations.txt" $CUR_DIR/data_files/${EXECUTABLES[${CUR_EXECUTABLE}]}_${THREADS_COUNT}_iterations.txt
	  mv "data_file_time.txt" $CUR_DIR/data_files/${EXECUTABLES[${CUR_EXECUTABLE}]}_${THREADS_COUNT}_time.txt

	  ((THREADS_COUNT++))
	  DIMS_COUNT=10
    done
    ((CUR_EXECUTABLE++))
    THREADS_COUNT=1
    cd ${CUR_DIR}
done

    


