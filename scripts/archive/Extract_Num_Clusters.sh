echo "Alpha test:"
fns=$(ls ../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha*_beta0.5.txt)
for fn in $fns; 
	do 
	echo $fn
	grep "Num clusters" $fn; 
	echo 
	done

echo
echo "Beta test:"
fns=$(ls ../output/StringDB/string_synapse_interactions_combined_score.clixo_alpha0.5_beta*.txt)
for fn in $fns; 
	do 
	echo $fn
	grep "Num clusters" $fn; 
	echo 
	done
