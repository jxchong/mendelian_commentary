set -eo pipefail


DATESTRING=$1;




perl parse_genemap2.pl --genemap2 raw_download_${DATESTRING}/genemap2.txt --out ${DATESTRING}.genemap2.parsed.txt
perl parse_omimtxtZ_count_NGS_year.pl --in raw_download_${DATESTRING}/omim.txt.gz --out ${DATESTRING}.omimtxtZ.parsed.NGS.year.inheritance.txt
perl combine_omimtxtZ_genemap2_count_NGS_year.pl --genemap2 ${DATESTRING}.genemap2.parsed.txt --mim2gene raw_download_${DATESTRING}/mim2gene.txt --omimtxt ${DATESTRING}.omimtxtZ.parsed.NGS.year.inheritance.txt --out ${DATESTRING}.combinedOMIM.mentionsNGS.year.inheritance.txt

# Rscript --vanilla code/genemap2/commentary_graphs_OMIM.R ${DATESTRING} >> ${DATESTRING}.values.txt
