#for i in 1 2 3 
#	do 
#		for j in 4 5 6 7 8 9  
#		do 
#			nohup bash ./run_clixo.sh 0.$i 0.$j &
#		done 
#	done 
#
#for i in 4 5 6
#	do 
#		for j in 1 2 3 7 8 9 
#		do 
#			nohup bash ./run_clixo.sh 0.$i 0.$j &
#		done 
#	done 

for i in 6 7 8 9
	do 
		for j in 1 2 3 4 5 6 7 8 9 
		do 
			nohup bash ./run_clixo.sh 0.$i 0.$j &
		done 
	done 
