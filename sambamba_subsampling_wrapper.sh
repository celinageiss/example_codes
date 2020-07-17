#!/bin/bash

sambamba view -f bam -t 8 --subsampling-seed=0 -s $THRESHOLD $BAM_FILE -o $OUTPUT