if [ "$1" == 'jr' ];
	then 
	clixo_script="/gpfs/data01/glasslab/home/joreyna/projects/BNFO286/clixo_0.3/clixo"
elif [ "$1" == 'as' ];
	then
	clixo_script="/mnt/c/Users/Anubhav/Desktop/mhk7-clixo_0.3-0362bea/mhk7-clixo_0.3-0362bea/clixo"
fi 

alpha=$2
beta=$3
interactions_fn="../output/StringDB/string_synapse_interactions_combined_score.tsv"
clixo_fn="../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha${alpha}_beta${beta}.txt"
echo $clixo_fn
cmd="$clixo_script $interactions_fn $alpha $beta > $clixo_fn"
eval $cmd
