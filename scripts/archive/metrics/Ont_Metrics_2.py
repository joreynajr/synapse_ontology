
# for questions, contact karenmei@ucsd.edu

#Goal: To evaluate the custom hierarchy (Metric #2)

#Metric 2: To create a pipeline for testing for how well the custom ontology maps onto reference ontology using DDOT 
#How well does the model capture the known structure of the reference ontology? (evaluated by alignment to the Gene Ontology, i.e. how many GO terms significantly overlap gene modules in your model?)

#	   code for customizing ontologies: DDOT: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb

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


#Create dictionary of all the term_names and description (in this example: GO)
#i.e. example in this dictionary: key is 'GO:0005694'; value: 'chromosome'
#input: Term is a GO term; for example: 'GO:0005694'

def Find_GO_Term_Desc(Term):
	with open('goID_2_name.tab', 'r') as f:
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
def Generate_Ontology_File(Term, file_name):
	ndex_server ='http://public.ndexbio.org' 
	go_human = Ontology.from_ndex('http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf')
	ont= go_human.focus(Term)
	ont.to_table(file_name)
	return ont

#Convert the custom ontology into the desired descriptions and gene names (in this example: GO descriptions and gene IDs)
#the input to this function is ont, which is the output of Generate_Ontology_File
def Find_GO_Focus_GeneDict(ont):
	ont = ont.propagate(direction='forward', gene_term=True, term_term=False)
	terms_to_genes=ont.term_2_gene
	GO_Desc=[]
	gene_names=[]
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
	for key in syn_ont:
		pair=[]
		ont_gene_list=syn_ont[key]
		size_ont_genes=len(ont_gene_list)
		overlap_num=len(list(set(ont_gene_list)&set(test_gene_list)))
		pair.append(key)
		pair.append(overlap_num)
		pair.append(size_ont_genes)
		overlap.append(pair)
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

#Find the number of enriched modules
#function input: ont1, ont2
#ont1 (ontology 1); ont2 (ontology 2)
def Num_Enriched_Modules(ont1, ont2):
	chr_ont=Find_GO_Focus_GeneDict(ont1)
	test_ont=Find_GO_Focus_GeneDict(ont2)
	enriched_modules=[]
	for key in test_ont:
		test_gene_list=test_ont[key]
		if len(test_gene_list)!=0:
			enriched_terms=Find_Enrichment(ont1, chr_ont, test_gene_list)
			enriched_modules.append(enriched_terms)
		else:
			continue
	enriched_modules = [x for x in enriched_modules if x != []]
	enriched_modules = [item for sublist in enriched_modules for item in sublist]
	enriched_term_names=[]
	for item in enriched_modules:
		enriched_term_name=item[0]
		enriched_term_names.append(enriched_term_name)
	unique_enriched_term_names=set(enriched_term_names)
	#print (unique_enriched_term_names)
	print ('Num of Enriched Modules', len(unique_enriched_term_names))
	return 

#Generate custom ontology 1 from Chromatin branch from human GO
ontology_file=Generate_Ontology_File('GO:0000785', 'ont1.txt')
ont1=Ontology.from_table('ont1.txt')
#print (ont1)

#Generate custom ontology 2 from Chromosome branch from human GO
ont2_file=Generate_Ontology_File('GO:0005694', 'ont2.txt')
ont2=Ontology.from_table('ont2.txt')
#print (ont2)

Num_Enriched_Modules(ont1, ont2)

