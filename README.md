# Dcitri_genome_EVE_scripts_2021

This repository contains the python scripts used to detect endogenous viral elements (EVEs) in Diaphorina citri and can be applied to any arthropod genome of interest. These scripts are adapted from "ter Horst et al. 2019". Detailed methodology can be accessed in "Carlson et al. 2022". 

"Parse_EVE_XML.py" parses a BLASTx XML output (outfmt 5), removing duplicate and overlapping viral hits, retaining hits with higher bitscores. 

"Run_Drosophila_BLASTx.py" extracts filtered EVE hits from "Parse_EVE_XML.py" output and performs a reverse BLASTx search against the Drosophila melanogaster proteome (Uniprot proteome accession number UP000000803) and removes any EVEs with hits to the proteome. 

Final output should then be checked against NCBI nr protein database using BLASTp to filter hits manually and generate the final EVE list. 
