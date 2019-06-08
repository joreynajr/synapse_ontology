#for questions, contact karenmei@ucsd.edu

#Goal: To evaluate the custom hierarchy (Metric #3)

#How well does the model organize psychiatric disease genes? (your data-driven ontology will be used to provide gene sets for functional enrichment)
#output: the number of disease genes found in significantly enriched modules


#	   code for customizing ontologies: DDOT: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb

import sys 
import pandas as pd
import networkx as nx
import numpy as np
import os
import ddot
from ddot import Ontology
import csv
import scipy.stats as ss
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests

os.environ['KMP_DUPLICATE_LIB_OK']='True'

#Create dictionary of all GO_term_names and description
def Find_GO_Term_Desc(Term):
	#with open('goID_2_name.tab', 'r') as f:
	#print(mod_goID_fn)
	with open(mod_goID_fn, 'r') as f:
		GO_terms=[r for r in csv.reader(f, delimiter='\t')]

	term=[]
	desc=[]
	for i in range(len(GO_terms)):
		term.append(GO_terms[i][0])
		desc.append(GO_terms[i][1])

	dictionary=dict(list(zip(term, desc)))
	description=dictionary[Term]

	return description


#Generate custom ontology file
def Generate_Ontology_File(Term):
	ndex_server ='http://public.ndexbio.org' 
	go_human = Ontology.from_ndex('http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf')
	ont= go_human.focus(Term)
	ont.to_table('custom_ontology.txt')
	return ont


def Find_GO_Focus_GeneDict(ont):
	ont = ont.propagate(direction='forward', gene_term=True, term_term=False)
	terms_to_genes=ont.term_2_gene
	gene_names=[]
	GO_Desc = []
	for key in terms_to_genes:
		GO_name=Find_GO_Term_Desc(key)

		GO_Desc.append(GO_name)

		genes_in_terms=terms_to_genes[key]
		genes_list=[]
		for item in genes_in_terms:
			gene_name=ont.genes[item]
			genes_list.append(gene_name)
		gene_names.append(genes_list)
	
	GO_term_gene=list(zip(GO_Desc, gene_names))
	
	GO_dict=dict(GO_term_gene)
	return GO_dict



#Find the enriched terms within the custom ontology
#function input: ont (output of Generate_Ontology_File)
#               syn_ont (output of Find_GO_Focus_GeneDict)
#               test_gene_list (list of genes to be tested)
def Find_Enrichment(ont, syn_ont, test_gene_list):
	num_ont_genes=len(ont.genes)
	num_intersect_test_ont=len(set(ont.genes)&set(test_gene_list))
	overlap=[]

	# Getting the intersection between each ontology term and the test terms 
	for key in syn_ont:

		pair=[]

		ont_gene_list=syn_ont[key]
		size_ont_genes=len(ont_gene_list)
		overlap_num=len(list(set(ont_gene_list)&set(test_gene_list)))

		pair.append(key)
		pair.append(overlap_num)
		pair.append(size_ont_genes)
		overlap.append(pair)

	# Performing the hypergeometric test 
	p_values=[]
	for item in overlap:
		intersect_termgenes_testgenes=item[1]
		num_term_genes=item[2]
		p=hypergeom.sf(intersect_termgenes_testgenes-1, num_ont_genes,num_intersect_test_ont, num_term_genes)
		p_values.append(p)
	p_adj=multipletests(p_values,method='fdr_bh')


	boolean_p=p_adj[0]
	idx_true=[i for i, e in enumerate(boolean_p) if e != False]
	p_adj=p_adj[1]
	min_p=min(p_adj)
	idx_min_p=[i for i, e in enumerate(p_adj) if e == min_p]
	most_enriched_terms=[]
	for item in idx_min_p:
		term_padj_pair=[]
		enriched=overlap[item][0]
		GO_term_padj=p_adj[item]
		term_padj_pair.append(enriched)
		term_padj_pair.append(GO_term_padj)
		most_enriched_terms.append(term_padj_pair)
	true_terms=[]

	for item in idx_true:
		true_pair=[]
		true_term=overlap[item][0]
		true_padj=p_adj[item]
		true_pair.append(true_term)
		true_pair.append(true_padj)
		true_terms.append(true_pair)
	#print ('terms above threshold:', true_terms)
	return true_terms

def Find_num_genes_in_enriched(ont, chr_ont, test_gene_list):

	true_terms=Find_Enrichment(ont, chr_ont, test_gene_list)
	genes_in_true_terms=[]

	for item in true_terms:
		term=item[0]
		genes=chr_ont[term]
		genes_in_true_terms.append(genes)
	print ('genes in true terms:', genes_in_true_terms)
	overlap=list(set(genes_in_true_terms[0])&set(test_gene_list))
	overlap_num=len(overlap)
	print (overlap_num)
	return overlap_num




## Test genes: a list of some genes involved in transcriptionally active chromatin
ont_fn = sys.argv[1]
genes_fn = sys.argv[2]

def replace_edgetype(e):
	if e == 'gene':
		return('Gene-Term')
	else:
		return('Child-Parent')

# Clean the ontology
mod_ont_fn = os.path.basename(ont_fn) + '.mod'
#with open(ont_fn) as f, open(mod_ont_fn, 'w') as fw: 
#
#	for line in f:
#		if line.startswith('#'):
#			continue 
#
#	header = line.split()[0:3]
#	fw.write('\t'.join(header) + '\n')	
#	for line in f: 
#		line = line.strip().split('\t')
#		line[2] = replace_edgetype(line[2])
#		record = '\t'.join(line[0:3])
#		fw.write(record + '\n')	
#clixo_ont = pd.read_table(ont_fn, dtype=str, comment='#', header=None, 
#                          names=['Parent', 'Child', 'EdgeType', 'drop'])
clixo_ont = pd.read_table(ont_fn, dtype=str, comment='#', header=0, 
                          names=['Parent', 'Child', 'EdgeType', 'drop'])
if clixo_ont.shape[0] == 4:
	clixo_ont.drop('drop', axis=1, inplace=True)
clixo_ont.loc[:, 'EdgeType'] = clixo_ont.EdgeType.str.replace('gene', 'Gene-Term').replace('default', 'Child-Parent')
clixo_ont.to_csv(mod_ont_fn, sep='\t', index=None)

# Make a dummy goID file  
mod_goID_fn = os.path.basename(ont_fn) + '.goID'
with open(mod_goID_fn, 'w') as fw: 
	uniq_parents = sorted(clixo_ont['Parent'].unique())
	for parent in uniq_parents: 
		fw.write('{}\t{}\n'.format(parent, parent))
		
#
#	for line in f:
#		if line.startswith('#'):
#			continue 
#
#	header = line.split()[0:3]
#	parents = set() 
#	for line in f: 
#		line = line.strip().split('\t')
#		parents.add(line[0])	
#
#	for parent in sorted(parents):
#		fw.write('{}\t{}'.format(parent, parent))

with open(genes_fn) as f: 
	test_gene_list = [x.strip() for x in f.readlines()] 

ont = Ontology.from_table(mod_ont_fn)

# Make dictionary of terms and genes within custom ontology
chr_ont = Find_GO_Focus_GeneDict(ont)

#print(chr_ont)

# Find number of test genes in enriched modules:
Find_num_genes_in_enriched(ont, chr_ont, test_gene_list)











