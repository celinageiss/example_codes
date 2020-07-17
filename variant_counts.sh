#!/bin/bash

SHORT_PID='$1'  # List of short patient IDs (without prefix)
COVERAGES='30x 60x 90x max'

OUTPUT_FILE='$2'  # path to output csv-file

# Header of output file
echo -e '#PID_group\tPID\ttumor_cell_content\tfunctional_SNVs\tfunctional_indels\tstructural_variants\tCN_ratios' > $OUTPUT_FILE

for PID in $SHORT_PID; do
	echo -e 'T021-'$PID >> $OUTPUT_FILE  # PID group
	
	# WGS files with different coverages
	for COV in $COVERAGES; do
	
		TCC=$(tail -n4 /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-${COV}/mpileup/snvs_T021-${PID}-${COV}_purityEST.txt | head -n1)

		SNVS=$(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-${COV}/mpileup/snvs_T021-${PID}-${COV}_somatic_functional_snvs_conf_8_to_10.vcf | wc -l)
		
		INDELS=$(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-${COV}/platypus_indel/indel_T021-${PID}-${COV}_somatic_functional_indels_conf_8_to_10.vcf | wc -l)
		
		SVS=$(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-${COV}/SOPHIA/svs_T021-${PID}-${COV}_filtered_somatic_minEventScore5.bedpe | wc -l)

		CNRATIOS=$(ls /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-${COV}/ACEseq/*ALL.png | cut -d'/' -f13 | cut -d'_' -f3) # Extracts possible ratios
		
		# write genome counts in output file
		echo -e '\tT021-'${PID}-${COV}'\t'$TCC'\t'$SNVS'\t'$INDELS'\t'$SVS'\t'$CNRATIOS >> $OUTPUT_FILE
	
	done
	
	# WES files
	
	TCC=$(tail -n4 /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-exome/mpileup/snvs_T021-${PID}-exome_purityEST.txt | head -n1)
	
	SNVS=$(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-exome/mpileup/snvs_T021-${PID}-exome_somatic_functional_snvs_conf_8_to_10.vcf | wc -l)
	
	INDELS=$(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-exome/platypus_indel/indel_T021-${PID}-exome_somatic_functional_indels_conf_8_to_10.vcf | wc -l)
	
	SVS=$(($(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-exome/crest/tumor_T021-${PID}-exome.DELDUPINV | grep 'somatic' | wc -l) + $(grep -v ^# /icgc/dkfzlsdf/analysis/hipo/hipo_021/illumina_comparison/downsampled_files/analysis/results_per_pid/T021-${PID}-exome/crest/tumor_T021-${PID}-exome.TX | grep 'somatic' | wc -l)))
	
	# write exome counts in output file
	echo -e '\tT021-'${PID}-exome'\t'$TCC'\t'$SNVS'\t'$INDELS'\t'$SVS >> $OUTPUT_FILE
	
done