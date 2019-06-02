clixo_script="/gpfs/data01/glasslab/home/joreyna/projects/BNFO286/clixo_0.3/clixo"
alpha=$1
beta=$2
interactions_fn="../output/StringDB/string_synapse_interactions_combined_score.tsv"
clixo_fn="../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha${alpha}_beta${beta}.txt"
echo $clixo_fn
cmd="$clixo_script $interactions_fn $alpha $beta > $clixo_fn"
eval $cmd
