#!/usr/bin/env perl
#
# Description:
#
#
#
# Created by Jessica Chong on 2015-02-22.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($morbidmapfile, $genemap2file, $outputfile, $help);

GetOptions(
	# 'morbidmap=s' => \$morbidmapfile,
	'genemap2=s' => \$genemap2file,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

# if (!defined $morbidmapfile) {
# 	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --morbidmap not defined.\n")
# } els
if (!defined $genemap2file) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --genemap2 not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
}




# For the file genemap2, the fields are, in order:
# 1 Chromosome (NCBI)
# 2 Genomic position start * (NCBI)
# 3 Genomic position end (NCBI)
# 4 Cyto location (OMIM)
# 5 Computed cyto location (UCSC)
# 6 MIM Number for Gene/Locus (OMIM)
# 7 Gene symbols (OMIM)
# 8 Gene name (OMIM)
# 9 Approved gene symbol (HGNC)
# 10 Entrez gene ID (NCBI)
# 11 Ensembl gene ID (Ensembl)
# 12 Comments (OMIM)
# 13 Phenotype(s) (OMIM)
# 14 Mouse gene symbol & ID (MGI)

my %genemap2data;
open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "phenoname\tphenoMIMnum\tphenoMappingKey\tLocusSymbols\tGeneMIMnum\tCytoLoc\tisComplex\tinheritDN\tinheritAD\tinheritAR\tinheritX\tinheritY\tinheritMT\tinheritMosaic\n";
open (my $input_handle, "$genemap2file") or die "Cannot read $genemap2file: $!.\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
    if ($_ =~ /^#/) {
        next;
    }
	my @values = split(/\t/, $_);
	if (!defined $values[12]) {		# if entry is gene only
		next;
	}
	my @splitphenos = split("; ", $values[12]);
	my ($chromosome, $locussymbol, $locusMIM, $cytoloc) = ($values[0], $values[8], $values[5], $values[3]);
	# my ($phenoname, $locussymbol, $locusMIM, $cytoloc) = split(/\t/, $_);

	foreach my $phenoname (@splitphenos) {
		my ($phenoMIM, $phenomappingkey) = qw(NA -9);
		my ($inheritance, $inheritDN, $inheritAD, $inheritAR, $inheritX, $inheritY, $inheritMT, $inheritMosaic) = (0) x 8;
		if ($phenoname =~ m/(\d{6})/) {
			$phenoMIM = $1;
		}
		if ($phenoname =~ m/\(([1-4])\)/) {
			$phenomappingkey = $1;
		}
		# if ($phenoname =~ m/\d{6} \(\d\), (.+)$/) {			# technically this is where the model should be in the entry name but it isn't always. I think it's probably fine to just search the entire phenotype name
			# $inheritance = $1;
			if ($phenoname =~ m/autosomal dominant/i) {
				$inheritAD = 1;
			}
			if ($phenoname =~ m/autosomal recessive/i) {
				$inheritAR = 1;
			}
			if ($phenoname =~ m/X-linked/i) {
				$inheritX = 1;
			}
			if ($phenoname =~ m/Y-linked/i) {
				$inheritY = 1;
			}
			# genemap2 doesn't contain any mitochondrial genes!
			# if ($phenoname =~ m/Mitochondrial/i && $chromosome eq 'MT') {
			# 	$inheritMT = 1;
			# }
			if ($phenoname =~ m/somatic|somatic mutation|mosaic[ism]* for|germline mosaic/i) {
				$inheritMosaic = 1;
			}
		# }
		if ($phenomappingkey == 4) {
			# if a chromosomal del/dup syndrome, the gene MIM *is usually* the phenotype MIM
			if ($phenoMIM eq 'NA')  {
				$phenoMIM = $locusMIM;
			}
		}
		# if no phenotype MIM number, the gene MIM *is usually* the phenotype MIM
		# if ($phenoMIM eq 'NA')  {
		# 	$phenoMIM = $locusMIM;
		# }
		my $isComplex = "no";
		# remove the following phenotypes:
		# QTL or quantitative trait locus
		# suscep* (susceptibility but susceptibility is spelled incorrectly in OMIM in a few entries)
		# risk
		# [] or {}
		# !!! do not remove somatic but should flag
		if ($phenoname =~ m/risk|QUANTITATIVE TRAIT LOCUS|QTL|multifactorial|suscep(\w+) to/i || $phenoname =~ /\[/ || $phenoname =~ /\{/) {
			$isComplex = "yes";
		}
		if ($phenoname =~ "somatic" && $phenoname =~ /carcinoma|cancer|tumor|leukemia|lymphoma|sarcoma|blastoma|adenoma|cytoma|myelodysplastic|Myelofibrosis|oma,/i ) {
			$isComplex = "cancer";
		} elsif ($phenoname =~ "somatic" && $isComplex ne 'yes') {
			$isComplex = "somatic";
		}

		# $genemap2data{$phenoMIM} = "$phenoname\t$phenoMIM\t$phenomappingkey\t$locussymbol\t$locusMIM\t$cytoloc\t$isComplex\t$inheritDN\t$inheritAD\t$inheritAR\t$inheritX\t$inheritY\t$inheritMT\t$inheritMosaic\n";
		print $output_handle "$phenoname\t$phenoMIM\t$phenomappingkey\t$locussymbol\t$locusMIM\t$cytoloc\t$isComplex\t$inheritDN\t$inheritAD\t$inheritAR\t$inheritX\t$inheritY\t$inheritMT\t$inheritMosaic\n";
		# print STDERR "genemap2\t$phenoMIM\n";
	}
}
close $input_handle;



close $output_handle;






################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


parse_morbidmap.pl - Parse morbidmap into a more useful form!


=head1 SYNOPSIS


perl B<xxxx.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--genemap2> F<input file>

	input file

=item B<--out> F<output file>

	name of output file

=item B<--help> I<help>

	print documentation

=back


=head1 FILES


xx


=head1 EXAMPLES


xxxxxxxx


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
