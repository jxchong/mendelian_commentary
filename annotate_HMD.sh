paste <(head -1 HMD_HumanPhenotype.rpt.header.txt) <(echo -e "MortalityBool\tNormalBool\tNoAnalysisBool\tHPOcount") > HMD_HumanPhenotype.rpt.mortality_HPOcount.txt

# mortality/aging MP:0010768
# no phenotypic analysis MP:0003012
# normal phenotype MP:0002873

perl -anF"\t" -e '
$F[$#F] =~ s/\s+$//;
print join("\t", @F);
if ($F[6] =~ "MP:0010768") {print "1";} else {print "0";}
if ($F[6] =~ "MP:0002873") {print "\t1";} else {print "\t0";}
if ($F[6] =~ "MP:0003012") {print "\t1";} else {print "\t0";}
my @ids = split(" ", $F[6]);
print "\t".scalar(@ids)."\n";
' HMD_HumanPhenotype.rpt.txt >> HMD_HumanPhenotype.rpt.mortality_HPOcount.txt
