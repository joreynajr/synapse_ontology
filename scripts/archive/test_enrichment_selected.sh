fns='../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.1_beta0.7.txt'
fns=$fns' ../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.1_beta0.9.txt'
fns=$fns' ../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.2_beta0.7.txt'
diseases="schizophrenia mdd bipolar autism adhd"
autism_fn=$1

for fn in $fns;  
	do 
		res=$(python test_gene_list_with_ontology.py $fn $autism_fn)
		echo $fn $res >> test_enrichment.$disease.results 
	done
