# Goal: To evaluate the custom hierarchy: #Metric 1
#	   Metric #1: How well does the model capture novel synapse proteins? (some of which have been recently uncovered by our collaborators using AP/MS/MS)
#	   code for customizing ontologies: DDOT: https://github.com/michaelkyu/ddot/blob/master/examples/Tutorial.ipynb

# JR Note: Need a strategy which uses a subbet of synapse genes 

import sys
import csv

import numpy as np
import pandas as pd

import networkx as nx

import scipy.stats as ss
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests

from matplotlib import pyplot as plt
import matplotlib; #matplotlib.use("TKAgg")

from ddot import Ontology
from igraph import *

os.environ['KMP_DUPLICATE_LIB_OK']='True'

#-----------------Metric 1------------------------------------------------
#define jaccard similarity index

def jaccard(a, b):
	s1 = set(a)
	s2 = set(b)
	return len(s1.intersection(s2)) / len(s1.union(s2))


def metric_1(ont_file, test_gene_list):
	ont1 = Ontology.from_table(ont_file)
	ont1_genes = ont1.genes
	test_recovery = jaccard(ont1_genes, test_gene_list)
	print ('recovery of test genes:', test_recovery)
	return test_recovery


# Example 1: toy ontology
# Toy example ontology
ont_file = '../data/toy_ontology.txt' 
# Gene list for toy ontology:
test_gene_list = sys.argv[2].split(',')
print (metric_1(ont_file, test_gene_list))


# Example 2: GO chromatin
def Generate_Ontology_File(term, file_name, url='http://public.ndexbio.org/v2/network/16565851-4bfc-11e9-9f06-0ac135e8bacf'):
	"""
	Generate custom ontology file.

	"""
	ndex_server ='http://public.ndexbio.org' 
	go_human = Ontology.from_ndex(url)
	ont= go_human.focus(term)
	ont.to_table(file_name)
	return ont

#Generate custom ontology from Chromatin branch from human GO
chr_ont = Generate_Ontology_File('GO:0000785', 'ont1.txt')
ont_file = 'ont1.txt'

#test_gene_list for chromatin GO:
test_gene_list = ['ANKRD17', 'ANKRD2', 'ANP32E', 'APTX', 'AR', 'ARID1A', 'ARID1B', 'ARRB1']
print (metric_1(ont_file, test_gene_list))

# Commandline format  
#ont_file = sys.argv[1]
#test_gene_list = sys.argv[2].split(',')
#print (metric_1(ont_file, test_gene_list))
