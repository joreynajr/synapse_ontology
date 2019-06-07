#for questions, contact karenmei@ucsd.edu

#Goal: To evaluate the custom hierarchy (Metric #3)

#How well does the model organize psychiatric disease genes? (your data-driven ontology will be used to provide gene sets for functional enrichment)
#output: the number of disease genes found in significantly enriched modules


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

#Create dictionary of all GO_term_names and description
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
def Generate_Ontology_File(Term):
	ndex_server ='http://public.ndexbio.org' 
	go_human = Ontology.from_ndex('http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf')
	ont= go_human.focus(Term)
	ont.to_table('custom_ontology.txt')
	return ont


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


#Test genes: a list of some genes involved in transcriptionally active chromatin

test_gene_list=['AFF4', 'PSIP1', 'H2AFB1', 'H2AFB2', 'H2AFB3', 'TTC37', 'WDR61', 'KMT2E', 'BCAS3', 'HIST1H1C', 'ZC3H8', 'CTR9', 'PADI2', 'PAF1', 'SUPT6H', 'EXOSC4', 'ELL', 'PCID2', 'EXOSC5', 'PELP1', 'ESR1', 'EXOSC10', 'EXOSC3', 'ICE1', 'ICE2']

test_gene_list = 'AAK1 ABCC8 ABHD17A ABHD17B ABHD17C ABHD6 ABI1 ABI2 ABL1'.split()

#Generate custom ontology from Chromatin branch from human GO
#ontology_file=Generate_Ontology_File('GO:0000785')
ont=Ontology.from_table('./ont2.txt')

#Make dictionary of terms and genes within custom ontology
chr_ont=Find_GO_Focus_GeneDict(ont)

#Find number of test genes in enriched modules:
Find_num_genes_in_enriched(ont, chr_ont, test_gene_list)

