# Pre-existing datasets to downloads

1. Human genes with mouse phenotypes for ortholog - HMD_HumanPhenotype.rpt.mortality_HPOcount.txt
	* Download data from MGI: http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
	* Run script `./annotate_HMD.sh ` to count genes with no phenotypic analysis, only normal phenotype reported, or linked to mortality/aging (mostly lethal)

2. HGNC gene list - HGNC_genes.tsv
	* Follow instructions here: https://github.com/macarthur-lab/gene_lists/blob/master/src/create_universe.bash

3. Constrained coding regions in 90%ile or higher - ccrs.autosomes.90orhigher.v2.20180420.bed.gz and ccrs.xchrom.90orhigher.v2.20180420.bed.gz
	* Download https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.90orhigher.v2.20180420.bed.gz
	* Download https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.90orhigher.v2.20180420.bed.gz
	* See https://github.com/quinlan-lab/ccr and/or Havrilla JM, Pedersen BS, Layer RM, Quinlan AR. 2018. A map of constrained coding regions in the human genome. Nat Genet 75: 1–12 for reference/more info

4. Genes depleted for nonsense-mediated decay-escaping variants - NMDescape.txt
	* Download Table S4 from Coban Akdemir Z, White JJ, Song X, Jhangiani SN, Fatih JM, Gambin T, Bayram Y, Chinn IK, Karaca E, Punetha J, Poli C, Baylor-Hopkins Center for Mendelian Genomics, Boerwinkle E, Shaw CA, Orange JS, Gibbs RA, Lappalainen T, Lupski JR, Carvalho CMB. 2018. Identifying Genes Whose Mutant Transcripts Cause Dominant Disease Traits by Potential Gain-of-Function Alleles. Am J Hum Genet 103: 171–187. (Free version here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6081281/)
	* Extract/copy first 3 columns to text file NMDescape.txt
	* Manually correct the first gene name to "MARCH11" (the file available from AJHG was mangled by Excel)

5. Request access for the following files from OMIM. Substitute your personal download/access token into the following URLs in place of "XXXXXXXXX"
```
wget --no-check-certificate https://omim.org/static/omim/data/mim2gene.txt
wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/mimTitles.txt
wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/genemap.txt
wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/morbidmap.txt
wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/genemap2.txt
wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/omim.txt.gz
```
