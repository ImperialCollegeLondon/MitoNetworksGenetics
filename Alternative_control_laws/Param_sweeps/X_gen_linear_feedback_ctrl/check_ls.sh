#!/bin/bash
for i in `seq 0 120`;
do
	printf -v DATA_DIR "p%d" $i	
	echo $DATA_DIR
	ls $DATA_DIR | wc -l
done
