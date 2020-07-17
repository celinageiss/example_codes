# Code Examples

This repository contains a couple of examples for codes I wrote during some of my former projects.

* _CIViC_annotation.py_  
This script takes a tsv table with mutated genes of a cancer patient as input and annotates hits for the identified variants in the database [CIViC](https://civicdb.org/home) using their API. It was implemented into the cancer patient anaylis pipeline of the German Cancer Research Centre (DKFZ) to identify possible tailored therapy options.

* _OT-2_PCR_purification.py_  
This programme was written for an automation of T cell receptor cloning with the liquid handling robot [OT-2](https://opentrons.com/ot-2/) from opentrons using their [API](https://docs.opentrons.com/v2/). It takes a 96-well plate and performs a bead-based PCR purification of the samples.

* _sambamba_subsampling_wrapper.sh_ and _sambamba_subsambpling_run.sh_  
These two scripts were used to perform a coverage downsampling of paired tumor-control whole-genome sequencing (WGS) data. The aim of the overall project was to compare the performance of deep whole-exome sequencing and different depths of WGS to find the best cost-effectiveness tradeoff for cancer diagnostics. The wrapper takes a table with patient IDs and parameters as an input and starts a cluster session where the sambamba-run script is called for each sample.

* _variant_comparison.sh_  
This script was used in the same project as the aforementioned sambamba scripts. It creates a summary output table for all samples, including the number of called variants (SNVs, indels, structural variants) for the different coverages of WGS and WES data.
