#!/bin/bash

# Downsampling of control
 
COV_TABLE="$1"  # contains coverages and TCC
PID="$2"  # input bam that needs to be downsampled
TUMOR_OR_CTRL="$3"  # t/T or c/C

DOWNSAMPLE_DIR="$4"  # output base directory

# Read COV_TABLE and create folder structure

grep -P "^$PID\t" $COV_TABLE | while read LINE
do
	TARGET_COV=$(echo "$LINE" | cut -f6)
	mkdir -p $DOWNSAMPLE_DIR/T021-${PID}-${TARGET_COV}x/alignment
	BAM_DIR=$DOWNSAMPLE_DIR/T021-${PID}-${TARGET_COV}x/alignment
	
	# Check if input is tumor or control

	if [ $TUMOR_OR_CTRL == "t" ] || [ $TUMOR_OR_CTRL == "T" ]
	then
		THRESHOLD=$(echo "$LINE" | cut -f8)
		BAM_FILE=$(echo "$LINE" | cut -f3)
		OUTPUT="$BAM_DIR/tumor_T021-${PID}-${TARGET_COV}x_rmdup.bam"
		LOG_FILE="$BAM_DIR/tumor_T021-${PID}-${TARGET_COV}x.log"
		
	elif [ $TUMOR_OR_CTRL == "c" ] || [ $TUMOR_OR_CTRL == "C" ]
	then
		THRESHOLD=$(echo "$LINE" | cut -f7)
		BAM_FILE=$(echo "$LINE" | cut -f2)
		OUTPUT="$BAM_DIR/control_T021-${PID}-${TARGET_COV}x_rmdup.bam"
		LOG_FILE="$BAM_DIR/control_T021-${PID}-${TARGET_COV}x.log"
		
	fi

	# request session on cluster and call wrapper_sambamba_subsampling
	qsub -o "$LOG_FILE" -j oe -l walltime=10:00:00,vmem=1G,mem=1G,nodes=1:ppn=8 -v "THRESHOLD=$THRESHOLD,BAM_FILE=$BAM_FILE,OUTPUT=$OUTPUT" /abi/data/illumina_comp/wrapper_sambamba_subsampling.sh
done