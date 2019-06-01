#for questions: karenmei@ucsd.edu

#Goal: To evaluate the custom hierarchy: #Metric 1
#	   Metric #1: How well does the model capture novel synapse proteins? (some of which have been recently uncovered by our collaborators using AP/MS/MS)

#	   code for customizing ontologies: DDOT: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb



import numpy as np
from igraph import *
import pandas as pd
import sys

sys.path.append("/Users/karenmei/Documents/ddot/")

import ddot
from ddot import Ontology
import matplotlib
matplotlib.use("TKAgg")
#print(matplotlib.get_backend())
from matplotlib import pyplot as plt
import networkx as nx
import csv
import scipy.stats as ss
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests

os.environ['KMP_DUPLICATE_LIB_OK']='True'

#-----------------Metric 1------------------------------------------------
#define jaccard similarity index

def jaccard(a, b):
	s1=set(a)
	s2=set(b)
	return len(s1.intersection(s2)) / len(s1.union(s2))


def metric_1(ont_file, test_gene_list):
	ont1=Ontology.from_table(ont_file)
	ont1_genes=ont1.genes
	test_recovery=jaccard(ont1_genes, test_gene_list)
	print ('recovery of test genes:', test_recovery)
	return test_recovery


#example 1: (toy ontology)---------------------------------------------------


#toy example ontology file:
ont_file='toy_ontology.txt'

#test_gene_list for toy ontology:
test_gene_list=['A', 'B', 'C', 'D', 'E']


print (metric_1(ont_file, test_gene_list))


#example 2: (GO chromatin)-------------------------------------------------

#Generate custom ontology file
def Generate_Ontology_File(Term, file_name):
	ndex_server ='http://public.ndexbio.org' 
	go_human = Ontology.from_ndex('http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf')
	ont= go_human.focus(Term)
	ont.to_table(file_name)
	return ont

#Generate custom ontology from Chromatin branch from human GO
chr_ont=Generate_Ontology_File('GO:0000785', 'ont1.txt')
ont_file='ont1.txt'

#test_gene_list for chromatin GO:
test_gene_list=['ANKRD17', 'ANKRD2', 'ANP32E', 'APTX', 'AR', 'ARID1A', 'ARID1B', 'ARRB1']


print (metric_1(ont_file, test_gene_list))

