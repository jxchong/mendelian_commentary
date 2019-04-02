This repository contains the scripts used to generate data for the manuscript:

Bamshad MJ, Nickerson DA, Chong JX. Mendelian gene discovery: reports of its demise have been greatly exaggerated.

With the files in this repository, you can generate the figures and the datasets used for analyses in the manuscript.


Instructions:


1. Download and format pre-existing datasets - `prep_data.md`
2. Analyze OMIM data to estimate year of discovery of the gene-disease relationship, year of syndrome delineation, and mode of inheritance `analyze_OMIM.sh DATE` where `DATE` is YYYY-MM-DD. The scripts assume that the OMIM data files are stored in a folder named "raw_download_DATE"
3. Generate all figures except Figure 2 `commentary_graphs_OMIM.R`
4. Generate Figure 2 `count_MC_candidate_genes.R` (graph drawing section should run interactively)
