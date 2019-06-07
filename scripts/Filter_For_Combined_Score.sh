 awk 'BEGIN{OFS="	"} {print $1, $2, $15}' ../output/StringDB/string_synapse_interactions.tsv > ../output/StringDB/string_synapse_interactions_combined_score.tsv
