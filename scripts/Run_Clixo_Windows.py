import sys 
import os
import subprocess

alpha = sys.argv[1]
beta = sys.argv[2]

clixo_script="/mnt/c/Users/Anubhav/Desktop/mhk7-clixo_0.3-0362bea/mhk7-clixo_0.3-0362bea/clixo"
interactions_fn = "../output/StringDB/string_synapse_interactions_combined_score.tsv"

# Naming of clixo output (can change if you need)
clixo_fn = "../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha{}_beta{}.txt".format(alpha, beta)

cmd = "{} {} {} {} > {}".format(clixo_script, interactions_fn, alpha, beta, clixo_fn)
print(cmd)
subprocess.call(cmd, shell=True)
