set -eo pipefail


DATESTRING=$1;

# requires requesting access for the following files from OMIM and obtaining a download key (your key shoudl be substituted in for XXXXXXXXX)


# wget --no-check-certificate https://omim.org/static/omim/data/mim2gene.txt
# wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/mimTitles.txt
# wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/genemap.txt
# wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/morbidmap.txt
# wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/genemap2.txt
# wget --no-check-certificate https://data.omim.org/downloads/XXXXXXXXX/omim.txt.gz


perl code/genemap2/parse_genemap2.pl --genemap2 raw_download_${DATESTRING}/genemap2.txt --out ${DATESTRING}.genemap2.parsed.txt
perl code/genemap2/parse_omimtxtZ_count_NGS_year.pl --in raw_download_${DATESTRING}/omim.txt.gz --out ${DATESTRING}.omimtxtZ.parsed.NGS.year.inheritance.txt
perl code/genemap2/combine_omimtxtZ_genemap2_count_NGS_year.pl --genemap2 ${DATESTRING}.genemap2.parsed.txt --mim2gene raw_download_${DATESTRING}/mim2gene.txt --omimtxt ${DATESTRING}.omimtxtZ.parsed.NGS.year.inheritance.txt --out ${DATESTRING}.combinedOMIM.mentionsNGS.year.inheritance.txt

Rscript --vanilla code/genemap2/commentary_graphs_OMIM.R ${DATESTRING} >> ${DATESTRING}.values.txt


# download the following files:
#
#
# HGNC_genes.tsv
# NMDescape.txt
# ccrs.autosomes.90orhigher.v2.20180420.bed.gz
# ccrs.xchrom.90orhigher.v2.20180420.bed.gz

# Human genes with mouse phenotypes for ortholog
# download Human-Mouse Phenotype list from MGI: http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt
# run annotate_HMD.sh to generate HMD_HumanPhenotype.rpt.mortality_HPOcount.txt

# HGNC gene list, see https://github.com/macarthur-lab/gene_lists/blob/master/src/create_universe.bash


# run count_MC_candidate_genes.R interactively in R to generate commentary_graphs_OMIM
