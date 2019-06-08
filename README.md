Making a synapse gene ontology.

StringDB
---
To collect the synapse network from StringDB we...After that we gather a list of 
disease genes from the OMIM database + those genes given to us and extract the 
first neighbors. 
1) Run the Extracting_String_Interaction.ipynb

Running the CliXo Pipeline
---
Once the network was build we ran clixo <network file> <alpha> <beta>

2) Open and execute Run_Clixo.ipynb (update Clixo path if necessary). 
 

Analysis
---
To be understand out ontology we performed a few analysis looking for: 
1) Visualize the new ontology.
2) Recall with respect to the GO ontology 
3) Enrichment of disease genes 
4) Analyze the number of children per parent node using MakeHistogram.ipynb 
