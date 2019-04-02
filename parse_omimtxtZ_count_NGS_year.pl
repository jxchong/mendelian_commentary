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


my $debuglevel = 0;

my ($inputfile, $outputfile, $help);

GetOptions(
	'in=s' => \$inputfile,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $inputfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
}



open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
# print $output_handle "MIMnum\tmimtitle\tMIMlink\tisComplex\tmentionsLinkage\tofficial_inheritDN\tofficial_inheritAD\tofficial_inheritAR\tofficial_inheritX\tofficial_inheritY\tofficial_inheritMT\tofficial_inheritMosaic\tofficial_inheritImprint\tinheritDN\tinheritAD\tinheritAR\tinheritX\tinheritY\tinheritMT\tinheritMosaic\tinheritImprint\tmentionsNGS\tmentionsNGSparagraph\tNGSyear\tyearDiscovered\tyearDelineated\n";
print $output_handle "MIMnum\tmimtitle\tMIMlink\tisComplex\tmentionsLinkage\t";
print $output_handle "official_inheritDN\tofficial_inheritAD\tofficial_inheritAR\tofficial_inheritX\tofficial_inheritY\tofficial_inheritMT\tofficial_inheritMosaic\tofficial_inheritImprint\t";
print $output_handle "inheritDN\tinheritAD\tinheritAR\tinheritX\tinheritY\tinheritMT\tinheritMosaic\tinheritImprint\t";
print $output_handle "mentionsNGS\tmentionsNGSparagraph\tNGSyear\tyearDiscovered\tyearDelineated\n";
my ($mentionsNGS, $NGSparagraph, $NGSyear, $yeardiscovered, $yeardelineated) = ('NA') x 5;
my ($official_inheritDN, $official_inheritAD, $official_inheritAR, $official_inheritX, $official_inheritY, $official_inheritMT, $official_inheritMosaic, $official_inheritImprint) = (0) x 8;
my ($mentionslinkage, $inheritDN, $inheritAD, $inheritAR, $inheritX, $inheritY, $inheritMT, $inheritMosaic, $inheritImprint) = (0) x 9;
my ($inrecord, $inTX, $inMISC, $inMG, $inCF, $inMAPPING, $inINHERITANCE, $inCSinherit) = (0) x 8;
my $printentry = '';
my (@MGtxt, @CFtxt, @MAPtxt, @INHERITtxt, @MISCtxt, @CS_Itxt, @TXtxt);
my ($MGtxt_paragraph, $CFtxt_paragraph, $MAPtxt_paragraph, $INHERITtxt_paragraph, $MISCtxt_paragraph, $CS_Itxt_paragraph, $TXtxt_paragraph) = ('') x 7;
my ($mimnum, $mimtitle);

my $input_handle;
# if ($inputfile =~ /\.Z$|\.gz$/) {
    # open ($input_handle, "zcat $inputfile |") or die "Cannot read $inputfile: $!.\n";
# } else {
    # open ($input_handle, "$inputfile") or die "Cannot read $inputfile: $!.\n";
# }


open ($input_handle, "zless $inputfile |") or die "Cannot read $inputfile: $!.\n";


while ( my $readline = <$input_handle> ) {
    if ($debuglevel >= 10) { print STDERR "readline = $readline"; }
	if (($readline =~ /^\*RECORD\*/ || $readline =~ /^\*THEEND\*/) && $inrecord == 1) {
		if ($debuglevel >= 1) { print STDERR "Begin new record\n"; }
		# if we are finishing a record/OMIM entry, then print stored info
		print $output_handle "$printentry\t";
		if ($mentionslinkage >= 1) {
			print $output_handle "1\t";
		} else {
			print $output_handle "0\t";
		}
		print $output_handle "$official_inheritDN\t$official_inheritAD\t$official_inheritAR\t$official_inheritX\t$official_inheritY\t$official_inheritMT\t$official_inheritMosaic\t$official_inheritImprint\t";
		print $output_handle "$inheritDN\t$inheritAD\t$inheritAR\t$inheritX\t$inheritY\t$inheritMT\t$inheritMosaic\t$inheritImprint\t";
		print $output_handle "$mentionsNGS\t$NGSparagraph\t$NGSyear\t$yeardiscovered\t$yeardelineated\n";

		($printentry, $mimnum, $mimtitle) = ('') x 3;
		$mentionsNGS = 'no';
		$NGSparagraph = $NGSyear = $yeardiscovered = 'NA';
		($MGtxt_paragraph, $CFtxt_paragraph, $MAPtxt_paragraph, $INHERITtxt_paragraph, $MISCtxt_paragraph, $TXtxt_paragraph, $CS_Itxt_paragraph) = ('') x 7;
		@MGtxt = @CFtxt = @MAPtxt = @INHERITtxt = @MISCtxt = @TXtxt = @CS_Itxt = ();
		$inTX = $inMISC = $inMG = $inCF = $inMAPPING = $inINHERITANCE = $inCSinherit = $inrecord = 0;
		($mentionslinkage, $inheritDN, $inheritAD, $inheritAR, $inheritX, $inheritY, $inheritMT, $inheritMosaic, $inheritImprint) = (0) x 9;
		($official_inheritDN, $official_inheritAD, $official_inheritAR, $official_inheritX, $official_inheritY, $official_inheritMT, $official_inheritMosaic, $official_inheritImprint) = (0) x 8;
	}

	if ($readline =~ /^\*RECORD\*/) {
		my $mimnum_expected = <$input_handle>;		# skip FIELD NO line
		my $mimnum = <$input_handle>;
		$mimnum =~ s/\s+$//;
		if ($mimnum_expected !~ /\*FIELD\* NO/) {
			die "line expected to contain *FIELD* NO, instead contains: $mimnum_expected\nnext line: $mimnum\n";
		}

        if ($debuglevel >= 1) { print STDERR "\n\nstarting entry #$mimnum\n"; }
		# if ($mimnum ne '118200') {
		#     last;
		# }
        # if ($mimnum != 100100) { next;}

		my $mimtitle_expected = <$input_handle>;		# skip FIELD TI line
		my $mimtitle = <$input_handle>;
		$mimtitle =~ s/\s+$//;
		if ($mimtitle_expected !~ /\*FIELD\* TI/) {
			die "line expected to contain *FIELD* TI, instead contains: $mimtitle_expected\nnext line: $mimtitle\n";
		}

		# An asterisk (*) before an entry number indicates a gene.
		# A number symbol (#) before an entry number indicates that it is a descriptive entry, usually of a phenotype, and does not represent a unique locus. The reason for the use of the number symbol is given in the first paragraph of the entry. Discussion of any gene(s) related to the phenotype resides in another entry(ies) as described in the first paragraph.
		# A plus sign (+) before an entry number indicates that the entry contains the description of a gene of known sequence and a phenotype.
		# A percent sign (%) before an entry number indicates that the entry describes a confirmed mendelian phenotype or phenotypic locus for which the underlying molecular basis is not known.
		# No symbol before an entry number generally indicates a description of a phenotype for which the mendelian basis, although suspected, has not been clearly established or that the separateness of this phenotype from that in another entry is unclear.
		# A caret (^) before an entry number means the entry no longer exists because it was removed from the database or moved to another entry as indicated.
		# Brackets, "[ ]", indicate "nondiseases," mainly genetic variations that lead to apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
		# Braces, "{ }", indicate mutations that contribute to susceptibility to multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).
		if ($mimtitle =~ /^[\*\^]/) {					# only keep entries that begin with #,%,+,no symbol; exclude if beginning with *, +, or ^
			$inrecord = 0;
			next;
		} else {
			$inrecord = 1;
		}
		my $isComplex = "no";

		if ($mimtitle =~ m/risk|QUANTITATIVE TRAIT LOCUS|QTL|multifactorial|suscep(\w+) to/i || $mimtitle =~ /\[/ || $mimtitle =~ /\{/ ) {
			$isComplex = "yes";
		}
		if ($mimtitle =~ /somatic/i && $isComplex ne 'yes') {
			$isComplex = "somatic";
		}
		$printentry =  "$mimnum\t$mimtitle\thttp://omim.org/entry/$mimnum\t$isComplex";
	}

	if ($readline =~ /^\*FIELD\* TX/i && $inrecord == 1) {			# don't move this block elsehwere or it will mess with recognition of the appropriate sections
	    $inTX = 1;
	    next;
	}
	if ($inTX == 1) {
        # print "in TX $readline\n";
	    if ($readline =~ /^- /i) { next; }
        if ($readline =~ /^[A-Z,_\-;: \/]+$/) {				# if beginning a new section of this OMIM entry (all caps title)
            if ($debuglevel >= 4) { print STDERR "TRIGGER end TEXT/TX section!!!!\n"; }
	        if ($TXtxt_paragraph !~ /^\s*$/) {			# exclude if line consists only of paragraph/line break with no content or paragraph only contains a link to another entry
	            push(@TXtxt, $TXtxt_paragraph);
	        }
			my $TXtxt_string = join(" ", @TXtxt);
			checkInheritance($TXtxt_string, \$inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint);
            $inTX = 0;
            if ($debuglevel >= 4) { print STDERR "ending TX\ncurrent = $readline\n\n"; }
            if ($debuglevel >= 4) { print STDERR "\nTX = $mentionslinkage\t@TXtxt\n\n"; }
	    } elsif ($readline =~ /^\n/ && $TXtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
	        if ($debuglevel >= 4) { print STDERR "saving new paragraph: $readline into $TXtxt_paragraph!!!!!!!!!!!!!! \n"; }
	        push(@TXtxt, $TXtxt_paragraph);
	        $TXtxt_paragraph = "";
	    } else {
	        $readline =~ s/\v//;
	        $TXtxt_paragraph .= " $readline";
	        if ($debuglevel >= 6) { print STDERR "add line $readline to $TXtxt_paragraph!!!!!!!!!!!!!!\n"; }
	    }
	}

	if ($readline =~ /^CLINICAL FEATURES$/i && $inrecord == 1) {
		$inCF = 1;
		$inTX = 0;
		next;
	}
	if ($inCF == 1) {
		if ($readline =~ /^- /i) { next; }
		if ($readline =~ /^[A-Z,_\-;: \/]+$/) {				# if beginning a new section of this OMIM entry (all caps title)
            if ($debuglevel >= 2) { print STDERR "TRIGGER end CF!!!! at $readline\n";}
			if ($CFtxt_paragraph !~ /^\s*$/) {			# if paragraph/line break exists with no content or paragraph only contains a link to another entry, do not count as a paragraph
				push(@CFtxt, $CFtxt_paragraph);
			}

			if ($debuglevel >= 7) { print STDERR print "CF paragraphs=".join("####", @CFtxt)."\n"; }
			if ($CFtxt[0] =~ /^\s*For [ \w]*[discussion|review] of [ ,\w]+, see [ \w(),]+\.$/) {		# if paragraph only contains a link to another entry, do not count as a paragraph
				my $test = shift(@CFtxt);
				# print "deleting $test\n";
			}

			$yeardelineated = 9999;
			$mentionslinkage += containsLinkage(join(" ", @CFtxt));

			if (scalar(@CFtxt) >= 1) {
				my @yearmatches;
				foreach my $para (@CFtxt) {					# find the earliest year in the entire Clinical Features section; future improvement: if the CF section has no year, then try the mapping section, example: https://omim.org/entry/605375
					@yearmatches = $para =~ /\([a-z, \.;-]*(19\d{2}|2\d{3})(?!\d)(?:[,;] [a-z, \.;-]*)*\)/ig;
					if (scalar(@yearmatches) >= 1) {
						foreach my $match (@yearmatches) {
							if ($match < $yeardelineated) {
								$yeardelineated = $match;
							}
						}
					}
				}
			}
			$inCF = 0;
            if ($debuglevel >= 2) { print STDERR "CF = @CFtxt\n\n"; }
		} elsif ($readline =~ /^\n/ && $CFtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
			if ($debuglevel >= 7) { print STDERR "saving new paragraph: $readline into $CFtxt_paragraph!!!!!!!!!!!!!! \n"; }
			push(@CFtxt, $CFtxt_paragraph);
			$CFtxt_paragraph = "";
		} else {
			$readline =~ s/\v//;
			$CFtxt_paragraph .= " $readline";
			if ($debuglevel >= 7) { print STDERR "add line $readline to $CFtxt_paragraph!!!!!!!!!!!!!!\n"; }
		}
	}

    # if ($readline =~ /^INHERITANCE$/i && $inrecord == 1) {			# don't move this block elsehwere or it will mess with recognition of the appropriate sections
    #     $inINHERITANCE = 1;
    #     $inCF = 0;
    #     next;
    # }
    # if ($inINHERITANCE == 1) {
    #     if ($debuglevel >= 3) { print STDERR "in INHERITANCE $readline\n";}
    #     if ($readline =~ /^- /i) { next; }
    #     if ($readline =~ /^[A-Z,_\-;: \/]+$/ || $readline =~ /^\*FIELD\*/) {				# if beginning a new section of this OMIM entry (all caps title)
    #         if ($debuglevel >= 3) { print STDERR "TRIGGER end INHERITANCE!!!!\n"; }
    #         if ($INHERITtxt_paragraph !~ /^\s*$/) {			# exclude if line consists only of paragraph/line break with no content or paragraph only contains a link to another entry
    #             push(@INHERITtxt, $INHERITtxt_paragraph);
    #         }
	# 		my $INHERITtxt_string = join(" ", @INHERITtxt);
	# 		$mentionslinkage += containsLinkage($INHERITtxt_string);
	# 		# don't try to get inheritance from this section when slurped in and analyzed section-wide because sometimes OMIM includes paragraphs like "xxx et al suggested that XX is dominantly inherited" (but xxx et al were incorrect)
	# 		# checkInheritance($INHERITtxt_string, \$inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint);
	# 		# if ($debuglevel >= 3) { print STDERR "inheritDN\tinheritAD\tinheritAR\tinheritX\tinheritY\tinheritMT\tinheritMosaic\tinheritImprint\n$inheritDN\t$inheritAD\t$inheritAR\t$inheritX\t$inheritY\t$inheritMT\t$inheritMosaic\t$inheritImprint\n$official_inheritDN\t$official_inheritAD\t$official_inheritAR\t$official_inheritX\t$official_inheritY\t$official_inheritMT\t$official_inheritMosaic\t$official_inheritImprint\n"; }
	#
	# 		$inINHERITANCE = 0;
    #         if ($debuglevel >= 3) { print STDERR "ending INHERITANCE\ncurrent = $readline\n\n"; }
    #         if ($debuglevel >= 3) { print STDERR "\nINHERITANCE = \t@INHERITtxt\n\n"; }
    #     } elsif ($readline =~ /^\n/ && $INHERITtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
    #         if ($debuglevel >= 3) { print STDERR "saving new paragraph: $readline into $INHERITtxt_paragraph!!!!!!!!!!!!!! \n"; }
    #         push(@INHERITtxt, $INHERITtxt_paragraph);
    #         $INHERITtxt_paragraph = "";
    #     } else {
    #         $readline =~ s/\v//;
    #         $INHERITtxt_paragraph .= " $readline";
    #         if ($debuglevel >= 7) { print STDERR "add line $readline to $INHERITtxt_paragraph!!!!!!!!!!!!!!\n"; }
    #     }
    # }

    if ($readline =~ /^MAPPING$/i && $inrecord == 1) {			# don't move this block elsehwere or it will mess with recognition of the appropriate sections
	    $inMAPPING = 1;
        $inINHERITANCE = 0;
	    next;
	}
	if ($inMAPPING == 1) {
        # print "in MAPPING $readline\n";
	    if ($readline =~ /^- /i) { next; }
        if ($readline =~ /^[A-Z,_\-;: \/]+$/) {				# if beginning a new section of this OMIM entry (all caps title)
            if ($debuglevel >= 4) { print STDERR "TRIGGER end mapping!!!!\n"; }
	        if ($MAPtxt_paragraph !~ /^\s*$/) {			# exclude if line consists only of paragraph/line break with no content or paragraph only contains a link to another entry
	            push(@MAPtxt, $MAPtxt_paragraph);
	        }
			my $MAPtxt_string = join(" ", @MAPtxt);
	        $mentionslinkage += containsLinkage($MAPtxt_string);
			# checkInheritance($MAPtxt_string, \$inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint);
            $inMAPPING = 0;
            if ($debuglevel >= 4) { print STDERR "ending MAPPING\ncurrent = $readline\n\n"; }
            if ($debuglevel >= 4) { print STDERR "\nMAPPING = $mentionslinkage\t@MAPtxt\n\n"; }
	    } elsif ($readline =~ /^\n/ && $MAPtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
	        if ($debuglevel >= 4) { print STDERR "saving new paragraph: $readline into $MAPtxt_paragraph!!!!!!!!!!!!!! \n"; }
	        push(@MAPtxt, $MAPtxt_paragraph);
	        $MAPtxt_paragraph = "";
	    } else {
	        $readline =~ s/\v//;
	        $MAPtxt_paragraph .= " $readline";
	        if ($debuglevel >= 6) { print STDERR "add line $readline to $MAPtxt_paragraph!!!!!!!!!!!!!!\n"; }
	    }
	}

	if ($readline =~ /^MOLECULAR GENETICS$/i && $inrecord == 1) {			# don't move this block elsehwere or it will mess with recognition of the appropriate sections
		$inMG = 1;
        $inMAPPING = 0;
		next;
	}
	if ($inMG == 1) {
		# print "in MG $readline\n";
		if ($readline =~ /Associations Pending Confirmation/i) { next; }
		if ($readline =~ /^- /i) { next; }

		if ($readline =~ /^\*FIELD\*/ || $readline =~ /^[A-Z,_\-;: \/]+$/) {
			if ($debuglevel >= 5) { print STDERR "detecting end of MG\n"; }
			if ($MGtxt_paragraph !~ /^\s*$/) {			# exclude lines that are paragraph/line break  with no content or paragraph only contains a link to another entry
				push(@MGtxt, $MGtxt_paragraph);
			}

			if ($debuglevel >= 5) { print STDERR "MGtxt=".join("####", @MGtxt)."\n"; }
			if ($MGtxt[0] =~ /^\s*For [ \w]*[discussion|review] of [ ,\w]+, see [ \w(),]+\.$/) {		# if paragraph only contains a link to another entry, do not count as a paragraph
				my $test = shift(@MGtxt);
				# print "deleting $test\n";
			}

			my $MGtxt_string = join(" ", @MGtxt);
			$mentionslinkage += containsLinkage($MGtxt_string);
			checkInheritance($MGtxt_string, \$inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint);
			if ($debuglevel >= 3) { print STDERR "inheritDN\tinheritAD\tinheritAR\tinheritX\tinheritY\tinheritMT\tinheritMosaic\tinheritImprint\n$inheritDN\t$inheritAD\t$inheritAR\t$inheritX\t$inheritY\t$inheritMT\t$inheritMosaic\t$inheritImprint\n$official_inheritDN\t$official_inheritAD\t$official_inheritAR\t$official_inheritX\t$official_inheritY\t$official_inheritMT\t$official_inheritMosaic\t$official_inheritImprint\n"; }
			if ($debuglevel >= 5) { print STDERR "MG section @MGtxt\n"; }

			for (my $i=0; $i<=$#MGtxt; $i++) {
				if ($mentionsNGS eq 'no') {
					if ($debuglevel >= 5) { print STDERR "paragraph $i\n"; }
					if ($debuglevel >= 5) { print STDERR "$MGtxt[$i]\n"; }
					($mentionsNGS, $NGSyear) = containsNGS($MGtxt[$i]);
					if ($mentionsNGS eq 'yes') {
						$NGSparagraph = "'".($i+1).'/'.($#MGtxt+1);
					}
				}
			}
			$yeardiscovered = 0;
			if (scalar(@MGtxt) >= 1) {
				# my @yearmatches = $MGtxt[0] =~ /\([a-z, \.-;]*(\d{4})(?:[,;] [\w, \.;-]*)*\)/ig;
				my @yearmatches = $MGtxt[0] =~ /\([a-z, \.;-]*(19\d{2}|2\d{3})(?!\d)(?:[,;] [a-z, \.;-]*)*\)/ig;
				if (scalar(@yearmatches) >= 1) {			# if length 1, then only one year
					foreach my $match (@yearmatches) {
						if ($match > $yeardiscovered) {
							$yeardiscovered = $match;
						}
					}
				}
			}
			$inMG = 0;
            if ($debuglevel >= 5) { print STDERR "2: $printentry\tmentionsNGS = $mentionsNGS\tparagraph = $NGSparagraph\n"; }
            if ($debuglevel >= 5) { print STDERR "\nMGtxt linkage = $mentionslinkage\n\n"; }
		} elsif ($readline =~ /^\n/ && $MGtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
            if ($debuglevel >= 5) { print STDERR "saving new paragraph: $MGtxt_paragraph!!!!!!!!!!!!!! \n"; }
			push(@MGtxt, $MGtxt_paragraph);
			$MGtxt_paragraph = "";
		} else {
			$readline =~ s/\v//;
			$MGtxt_paragraph .= " $readline";
            if ($debuglevel >= 5) { print STDERR "add line: $MGtxt_paragraph!!!!!!!!!!!!!!\n"; }
		}
	}

	if ($readline =~ /^INHERITANCE\:$/i && $inrecord == 1) {			# don't move this block elsehwere or it will mess with recognition of the appropriate sections
		$inCSinherit = 1;
		$inMG = 0;
		next;
	}
	if ($inCSinherit == 1) {
		# print "see CS MISC $readline\n";
		if ($readline =~ /^- /i) { next; }
		if ($readline =~ /^\*FIELD\*/ || $readline =~ /^[A-Z,_\-;: \/]+$/) {				# if beginning a new section of this OMIM entry (all caps title)
			if ($debuglevel >= 2) { print STDERR "TRIGGER end CS INHERITANCE!!!!\n"; }
			if ($CS_Itxt_paragraph !~ /^\s*$/) {			# exclude if line consists only of paragraph/line break with no content or paragraph only contains a link to another entry
				push(@CS_Itxt, $CS_Itxt_paragraph);
			}
			my $MISCtxt_string = join(" ", @CS_Itxt);
			checkInheritance($MISCtxt_string, \$inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint);
			checkInheritance($MISCtxt_string, \$official_inheritDN, \$official_inheritAD, \$official_inheritAR, \$official_inheritX, \$official_inheritY, \$official_inheritMT, \$official_inheritMosaic, \$official_inheritImprint);
			if ($debuglevel >= 2) { print STDERR "CS INHERITANCE section @MISCtxt\ninheritDN = $inheritDN\n"; }
			$inCSinherit = 0;
		} elsif ($readline =~ /^\n/ && $CS_Itxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
			push(@CS_Itxt, $CS_Itxt_paragraph);
			$CS_Itxt_paragraph = "";
		} else {
			$readline =~ s/\v//;
			$CS_Itxt_paragraph .= " $readline";
		}
	}
	if ($readline =~ /^MISCELLANEOUS\:$/i && $inrecord == 1) {			# don't move this block elsehwere or it will mess with recognition of the appropriate sections
		$inMISC = 1;
		$inCSinherit = 0;
		next;
	}
	if ($inMISC == 1) {
		# print "see CS MISC $readline\n";
		if ($readline =~ /^- /i) { next; }
		if ($readline =~ /^\*FIELD\*/ || $readline =~ /^[A-Z,_\-;: \/]+$/) {				# if beginning a new section of this OMIM entry (all caps title)
			if ($debuglevel >= 2) { print STDERR "TRIGGER end MISC!!!!\n"; }
			if ($MISCtxt_paragraph !~ /^\s*$/) {			# exclude if line consists only of paragraph/line break with no content or paragraph only contains a link to another entry
				push(@MISCtxt, $MISCtxt_paragraph);
			}
			$official_inheritDN += join(" ", @MISCtxt) =~ /de novo/gi;		# search for "de novo" and count matches in this section
			# my $MISCtxt_string = join(" ", @MISCtxt);
			# checkInheritance($MISCtxt_string, \$official_inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint);
			if ($debuglevel >= 2) { print STDERR "MISC section @MISCtxt\ninheritDN = $inheritDN\n"; }
			$inMISC = 0;
		} elsif ($readline =~ /^\n/ && $MISCtxt_paragraph !~ /^\s*$/) {				# if only content is a newline character, then the next line is the beginning of a new paragraph
			push(@MISCtxt, $MISCtxt_paragraph);
			$MISCtxt_paragraph = "";
		} else {
			$readline =~ s/\v//;
			$MISCtxt_paragraph .= " $readline";
		}
	}
}
close $input_handle;
close $output_handle;






sub containsNGS {
	my $text = $_[0];
	my $NGS = 'no';
	my $NGSyear = 0;

	if ($text =~ /exome[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /genome[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /massively[\s\-]*parallel[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /next[\s\-]*generation[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /high[\s\-]*throughput[\s\-]*sequencing/i) {
		$NGS = 'yes';
	}
	if ($text =~ /exome[\s\-]*capture/i) {
		$NGS = 'yes';
	}
	if ($text =~ /whole[\s\-]*exome/i) {
		$NGS = 'yes';
	}

	if ($NGS eq 'yes') {
		while ($text =~ /\((\d{4})\)/g) {
			if ($1 > $NGSyear) {
				$NGSyear = $1;
			}
		}
		if ($NGSyear < 2009) {
			# must be a false positive detection of NGS, like the mention of high-throughput sequencing in entry 243800 for Johanson-Blizzard Syndrome
			$NGS = 'no';
			$NGSyear = 'NA';
		}
	} else {
		$NGSyear = 'NA';
	}

	return ($NGS, $NGSyear);
}

sub containsLinkage {
	my $text = $_[0];
    my $mimnum = $_[1];

	my $linkage = 0;

	if ($text =~ /(\w+[\- ]generation)/i) {
		if ($1 !~ /next[\s\-]*generation[\s\-]*sequenc/i) {
			$linkage += 1;
            if ($debuglevel >= 6) { print STDERR "$text\nmatch $1\n"; }
		}
	}
	if ($text =~ /(linkage analysis|linkage mapping|lod score|lod \(|point linkage|linkage study)/i) {
		$linkage += 1;
        if ($debuglevel >= 6) { print STDERR "$text\nmatch $1\n"; }
	}
    if ($text =~ /(autosomal (dominant|recessive) inheritance)/i) {
        $linkage += 1;
        if ($debuglevel >= 6) { print STDERR "$text\nmatch $1\n"; }
    }
	return ($linkage);
}

sub checkInheritance {
	my $text = $_[0];
	my @values = (0) x 8;
	# $inheritDN, \$inheritAD, \$inheritAR, \$inheritX, \$inheritY, \$inheritMT, \$inheritMosaic, \$inheritImprint
	${$_[1]} += $text =~ m/de novo/gi;
	${$_[2]} += $text =~ m/autosomal dominant|haploinsufficiency/gi;
	${$_[3]} += $text =~ m/autosomal recessive|homozygous (|mutation)|recessive mutations|compound heterozygo/gi;
	${$_[4]} += $text =~ m/X-linked/gi;
	${$_[5]} += $text =~ m/Y-linked/gi ;
	${$_[6]} += $text =~ m/Mitochondrial/gi;
	${$_[7]} += $text =~ m/somatic mutation|mosaic(ism|) for|germline mosaic/gi;
	${$_[8]} += $text =~ m/imprinted gene|imprinting/gi;

	if ($text =~ m/heterozygous mutation|heterozygosity for (a \w+) mutation/i && $text !~ m/compound heterozygo|X-linked|mitochondrial/i) {
		${$_[2]} += 1;
	}
	if ($text =~ m/homozygous mutation|homozygosity for (a \w+) mutation/i && $text !~ m/X-linked|mitochondrial/i) {
		${$_[3]} += 1;
	}
}



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head NAME


xxx.pl -


=head SYNOPSIS


perl B<xxxx.pl> I<[options]>


=head ARGUMENTS


=over 4

=item B<--in> F<input file>

	input file

=item B<--out> F<output file>

	name of output file

=item B<--help> I<help>

	print documentation

=back


=head FILES


xx


=head EXAMPLES


xxxxxxxx


=head AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
