autism_fn='../output/omim_psychiatric_disease_genes/$1.txt'
fns='../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.1_beta0.7.txt'
fns=$fns' ../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.1_beta0.9.txt'
fns=$fns' ../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.2_beta0.7.txt'
for fn in $(ls ../output/StringDB/string_synapse_interactions_combined_score.clixo_*);
	do 
		res=$(python test_gene_list_with_ontology.py $fn $autism_fn)
		echo $fn $res >> test_enrichment.$1.results 
	done
