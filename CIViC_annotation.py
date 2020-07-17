#!/usr/bin/env python

# =============================================================================
# Name:     CIViC annotation
# Author:   Celina Geiss <celina.geiss@dkfz-heidelberg.de>
# Version:  1.3
# =============================================================================

import argparse
import urllib.request
import requests
import re
import json
import time


# Parser

parser = argparse.ArgumentParser(description = 'This program compares variants with entries of CIViC and PubMed and writes the number of matches in the output file. Accepts tables for integrated mutations, SNVs and indels as tsv-file. Execute with Python 3 and install the module "requests" if necessary.')

parser.add_argument('INPUT_FILE', help = 'Input file path')
parser.add_argument('OUTPUT_FILE', help = 'Output file path')

inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('-s', '--snvs', action = 'store_true', help = 'SNV table')
inputtype.add_argument('-i', '--indels', action = 'store_true', help = 'indels table')

parser.add_argument('-p', '--pubmed',  action = 'store_true', help = 'Annotate number of PubMed search results for [gene + neoplasms]. (Warning: Can slow down process significantly!)')

args = parser.parse_args()


# Input file path
input_file = str(args.INPUT_FILE)


# Output file path
output_file = str(args.OUTPUT_FILE)


# CIViC API
civic_url = 'https://civicdb.org/api/'


### Define functions

# Open website and convert to string
def convert_url_to_string(url):
	open_url = urllib.request.urlopen(url)
	url_bytes = open_url.read()
	url_string = url_bytes.decode('utf-8')
	return url_string


# Convert JSON to dictionary
def convert_json(url):
	url_string = convert_url_to_string(url)
	url_json = json.loads(url_string)
	return url_json


# Search for alternative gene name (alias) in CIViC
def alias(gene):
	json_page_1 = convert_json(civic_url + 'genes?count=25&page=1')
	total_pages = json_page_1['_meta']['total_pages']
	total_count = json_page_1['_meta']['total_count']

	civic_name = None

	# Iterate over all pages
	for i in range(1, total_pages + 1):
		open_gene_list_page = urllib.request.urlopen(civic_url + 'genes?count=25&page=' + str(i))
		gene_list_page_bytes = open_gene_list_page.read()
		gene_list_page = gene_list_page_bytes.decode('utf-8')

		json_page = json.loads(gene_list_page)
		entries_per_page = len(json_page['records'])

		# Iterate over all entries on a page
		for j in range(0, entries_per_page):
			aliases = json_page['records'][j]['aliases']
			if gene in aliases:
				civic_name = json_page['records'][j]['name']
				break
			else:
				pass

		if civic_name != None:
			break

	return civic_name


# Annotate PubMed results
def pubmed(gene):
	# PubMed API
	pubmed_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%s[TIAB]+AND+neoplasms[MeSH]'%(gene)

	pubmed_xml = convert_url_to_string(pubmed_url)
	pubmed_count = re.search(r'<Count>(\d+)</Count>', pubmed_xml).group(1)
	gene_dict[gene].append(pubmed_count)

	time.sleep(1) # PubMed cannot handle more than 3 requests per second, so pause for 1s

	return(gene_dict[gene])


# SNV lookup in CIViC
def snv_lookup(civic_entry):

	civic_variants = len(civic_entry['variants']) # how many variants were found for gene in CIViC

	gene_dict[gene].append(civic_variants) # gene_dict[gene][1] = total number of variants

	truncations = 0
	snvs = set()
	hits = [] # set of exact hits (truncations or position matches)

	# Look for variants
	for i in range(0, civic_variants):
		variant_name = civic_entry['variants'][i]['name']

		truncation = variant_name.find('TRUNCAT')
		if truncation >= 0:
			truncations += 1

		# SNVs
		snv = re.search(r'([A-Z]\d+)[A-Z].*?', variant_name)  # look for SNV at exact position

		if snv:
			snvs.add(snv.group(1))

	gene_dict[gene].append(len(snvs)) # gene_dict[gene][2] = total number of snvs


	for pos in gene_dict[gene][0]: # iterate over positions of input gene

		if gene_dict[gene][0][pos][1] == True and truncations > 0:
			hits.append('truncating_variant')

		for item in gene_dict[gene][0][pos][0]:
			if item in snvs:
				hits.append(item)

		gene_dict[gene][0][pos].append(sorted(hits))


# Indel lookup in CIViC
def indel_lookup(civic_entry):

	civic_variants = len(civic_entry['variants'])

	gene_dict[gene].append(civic_variants) # total number of variants found for gene in CIViC

	truncations = 0
	deletions = 0
	insertions = 0
	frameshifts = 0
	deletion_positions = []
	hits = []

	# Interate over variants for one gene
	for i in range(0, civic_variants): # iterate over all variants of one gene in CIViC
		variant_name = civic_entry['variants'][i]['name']

		truncation = variant_name.find('TRUNCAT')
		if truncation >= 0:
			truncations += 1

		if variant_name.find('INS') >= 0 or variant_name.find('ins') >= 0: # find insertions
			insertions += 1

		if variant_name.find('FRAME') >= 0 or variant_name.find('fs') >= 0: # find FRAMESHIFT or FRAME SHIFT
			frameshifts += 1

		if variant_name.find('DEL') >= 0 or variant_name.find('del') >= 0: # find deletions
			deletions += 1
			del_position = re.findall(r'\d+', variant_name)

			if len(del_position) > 1: # more than one base deleted
				del_range = [int(del_position[0]), int(del_position[1])+1]
				deletion_positions.append(del_range)

			elif len(del_position) == 1:
				del_range = [int(del_position[0]), int(del_position[0])+1]
				deletion_positions.append(del_range)


	for pos in gene_dict[gene][0]: # positions of mutations in input

		if gene_dict[gene][0][pos][1] is True and truncations > 0: # truncating mutation in input
			hits.append('truncating_variant')

		if len(deletion_positions) > 0: # find exact position hits
			for mut in gene_dict[gene][0][pos][0]:
				input_pos = mut[1:]
				for del_range in deletion_positions:
					if int(input_pos) in range(del_range[0], del_range[1]):
						interval = '-'.join([str(del_range[0]),str(del_range[1])])
						hits.append(interval)

		gene_dict[gene][0][pos].append(sorted(hits))

	num_indels = deletions + insertions + frameshifts + truncations

	gene_dict[gene].append(num_indels)



# =============================================================================
# SNVs
# =============================================================================

if args.snvs:

	with open(input_file, 'r') as snv_table_in:

		head = snv_table_in.readline().rstrip()
		head_split = head.split('\t')

		pos_column = head_split.index('POS')
		gene_column = head_split.index('GENE')
		annovar_column = head_split.index('ANNOVAR_TRANSCRIPTS')
		annovar_function_column = head_split.index('ANNOVAR_FUNCTION')
		exonic_classification_column = head_split.index('EXONIC_CLASSIFICATION')

		if args.pubmed:
			output_header = head + '\tCIViC_variant_entries\tCIViC_SNVs\tCIViC_exact_hits\tCIViC_gene_alias\tPubMed_entries\n'
		else:
			output_header = head + '\tCIViC_variant_entries\tCIViC_SNVs\tCIViC_exact_hits\tCIViC_gene_alias\n'

		output_list = []

		gene_column = head.split('\t').index('GENE')

		input_gene_list = []
		gene_dict = dict()

		for line in snv_table_in.readlines():
			line_split = line.rstrip().split('\t')

			if len(line_split) < len(head_split): # manche Zeilen haben eine Spalte weniger
				line = line.rstrip() + '\t.'

			output_list.append(line.rstrip())

			gene_raw = line_split[gene_column]
			gene = (gene_raw.split('(')[0]).split(',')[0] # manche Gene haben Zusatz in () hinter dem Namen --> extrahiere Name
			input_gene_list.append(gene)

			pos = line_split[pos_column]


			## Suche Input-Positionen
			input_variants = set(re.findall(r':p.([A-Z]\d+)', line_split[annovar_column]))

			if line_split[annovar_function_column] == 'splicing' or line_split[exonic_classification_column] == 'stopgain':
				truncation_in_input = True
			else:
				truncation_in_input = False

			if gene not in gene_dict:
				gene_dict[gene] = [dict()]

			gene_dict[gene][0][pos] = [input_variants, truncation_in_input]

	genes = ','.join(input_gene_list)

	civic_entry_url = ''.join([civic_url, 'genes/%s?identifier_type=entrez_symbol'%(genes)])
	civic_entries = convert_json(civic_entry_url) # Dict with all entries of genes that are in CIViC

	genes_in_civic = set() # Set with all genes that were found in CIViC
	not_in_civic = set()

	for entry in range(0,len(civic_entries)):

		civic_entry = civic_entries[entry]
		gene = civic_entry['name']
		genes_in_civic.add(gene)

		snv_lookup(civic_entry) # Call function for lookup in CIViC

		gene_dict[gene].append(0) # append 0 for alias

	not_in_civic = set(input_gene_list) - genes_in_civic

	# Search for aliases for genes that could not be found in CIViC
	for gene in not_in_civic:

		alternative_gene_name = alias(gene)

		if alternative_gene_name != None:
			gene_url = 'https://civic.genome.wustl.edu/api/genes/%s?identifier_type=entrez_symbol'%(alternative_gene_name)
			civic_entry = convert_json(gene_url)

			snv_lookup(civic_entry) # Call function for lookup in CIViC
			print(gene, alternative_gene_name)

			if gene != alternative_gene_name:
				gene_dict[gene].append(alternative_gene_name)
			else:
				gene_dict[gene].append(0) # append 0 for alias

		else:
			for pos in gene_dict[gene][0]:
				gene_dict[gene][0][pos].append([])

			gene_dict[gene].extend(['no entry in CIViC', 0, 0]) # no entry, no SNVs, no alias


	# Write in output file
	with open(output_file, 'w') as snv_table_out:

		snv_table_out.write(output_header)

		i = 0

		for gene in input_gene_list:

			# PubMed
			if args.pubmed:
				gene_dict[gene] = pubmed(gene)
				pubmed_count = gene_dict[gene][-1]
			else:
				pubmed_count = ''

			civic_count = gene_dict[gene][1]
			civic_snvs = gene_dict[gene][2]
			civic_alias = gene_dict[gene][3]

			pos = output_list[i].split('\t')[pos_column]
			hits = gene_dict[gene][0][pos][2]

			if len(hits) > 0:
				civic_hits = ','.join(hits)
			else:
				civic_hits  = 0

			output_line = '\t'.join([output_list[i], str(civic_count), str(civic_snvs), str(civic_hits), str(civic_alias), str(pubmed_count)+ '\n'])
			snv_table_out.write(output_line)

			i += 1


# =============================================================================
# Indels
# =============================================================================

elif args.indels:

	with open(input_file, 'r') as indels_table_in:

		head = indels_table_in.readline().rstrip()
		head_split = head.split('\t')

		pos_column = head_split.index('POS')
		gene_column = head_split.index('GENE')
		annovar_column = head_split.index('ANNOVAR_TRANSCRIPTS')
		annovar_function_column = head_split.index('ANNOVAR_FUNCTION')
		exonic_classification_column = head_split.index('EXONIC_CLASSIFICATION')

		if args.pubmed:
			output_header = head + '\tCIViC_variant_entries\tCIViC_indels\tCIViC_exact_hits\tCIViC_gene_alias\tPubMed_entries\n'
		else:
			output_header = head + '\tCIViC_variant_entries\tCIViC_indels\tCIViC_exact_hits\tCIViC_gene_alias\n'

		output_list = []

		gene_column = head.split('\t').index('GENE')

		input_gene_list = []
		gene_dict = dict()

		for line in indels_table_in.readlines():
			line_split = line.rstrip().split('\t')

			if len(line_split) < len(head_split): # manche Zeilen haben eine Spalte weniger
				line = line.rstrip() + '\t.'

			output_list.append(line.rstrip())

			gene_raw = line_split[gene_column]
			gene = (gene_raw.split('(')[0]).split(',')[0] # manche Gene haben Zusatz in () hinter dem Namen --> extrahiere Name
			input_gene_list.append(gene)

			pos = line_split[pos_column]

			# Look for exact input position
			input_variants = set(re.findall(r':p.([A-Z]\d+)', line_split[annovar_column]))

			if line_split[annovar_function_column] == 'splicing' or line_split[exonic_classification_column] == 'stopgain':
				truncation_in_input = True
			else:
				truncation_in_input = False

			if gene not in gene_dict:
				gene_dict[gene] = [dict()]

			gene_dict[gene][0][pos] = [input_variants, truncation_in_input]


	genes = ','.join(input_gene_list)
	civic_entry_url = ''.join([civic_url, 'genes/%s?identifier_type=entrez_symbol'%(genes)])
	civic_entries = convert_json(civic_entry_url) # Dict with all entries of genes that are in CIViC

	genes_in_civic = set() # Set with all genes that were found in CIViC
	not_in_civic = set()


	for entry in range(0,len(civic_entries)):

		civic_entry = civic_entries[entry]
		gene = civic_entry['name']
		genes_in_civic.add(gene)

		indel_lookup(civic_entry) # Call function for lookup in CIViC
		gene_dict[gene].append(0) # append 0 for alias

	not_in_civic = set(input_gene_list) - genes_in_civic


	# Search for aliases for genes that could not be found in CIViC
	for gene in not_in_civic:

		alternative_gene_name = alias(gene)

		if alternative_gene_name != None:
			gene_url = 'https://civic.genome.wustl.edu/api/genes/%s?identifier_type=entrez_symbol'%(alternative_gene_name)
			civic_entry = convert_json(gene_url)

			indel_lookup(civic_entry) # Call function for lookup in CIViC

			if gene != alternative_gene_name:
				gene_dict[gene].append(alternative_gene_name)
			else:
				gene_dict[gene].append(0) # append 0 for alias

		else:
			for pos in gene_dict[gene][0]:
				gene_dict[gene][0][pos].append([])
			gene_dict[gene].extend(['no entry in CIViC', 0, 0]) # no entry, no indels, no alias


	# Write in output file
	with open(output_file, 'w') as indels_table_out:

		indels_table_out.write(output_header)

		i = 0

		for gene in input_gene_list:

			# PubMed
			if args.pubmed:
				gene_dict[gene] = pubmed(gene)
				pubmed_count = gene_dict[gene][-1]
			else:
				pubmed_count = ''

			civic_count = gene_dict[gene][1]
			civic_indels = gene_dict[gene][2]
			civic_alias = gene_dict[gene][3]

			pos = output_list[i].split('\t')[pos_column]
			hits = gene_dict[gene][0][pos][2]

			if len(hits) > 0:
				civic_hits = ','.join(hits)
			else:
				civic_hits  = 0

			output_line = '\t'.join([output_list[i], str(civic_count), str(civic_indels), str(civic_hits), str(civic_alias), str(pubmed_count)+ '\n'])
			indels_table_out.write(output_line)

			i += 1
