Making a synapse gene ontology.

OMIM
---
Download the OMIM database file. 
1) Run OMIM_Psychiatric_Diseases.ipynb

StringDB
---
To collect the synapse network from StringDB we...After that we gather a list of 
disease genes from the OMIM database + those genes given to us and extract the 
first neighbors. 
2) Run the Extracting_String_Interaction.ipynb

Running the CliXo Pipeline
---
Once the network was build we ran clixo <network file> <alpha> <beta>

3) Open and execute Run_Clixo.ipynb (update Clixo path if necessary). 

Analysis
---
To be understand out ontology we performed a few analysis looking for: 
1) Visualize the new ontology using Visualize_Ontology.ipynb. 
2) Recall with respect to the GO ontology 
3) Enrichment of disease genes 
4) Analyze the number of children per parent node using MakeHistogram.ipynb 
