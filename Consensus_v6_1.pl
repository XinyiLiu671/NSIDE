#!/usr/bin/perl -w
use strict;
use File::Basename;
use Array::Utils qw[unique];
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Consensus.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to generate a consensus set of
NSIDEs.

------- ------- ------- ------- ------- ------- -------

Iutput files:

1. NSIDE.txt, ValNSIDE.txt, ConNSIDE.txt or
   DicedNSIDE.txt files, including at least 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Output file:

1. ConNSIDE.txt: consensus NSIDE, including 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Consensus.pl -i NSIDE*.txt [Options]

NSIDE_Consensus.pl -i NSIDE*.txt -n name
                   -o output_folder_name 

------- ------- ------- ------- ------- ------- -------

Options:
 -i|--input   > NSIDE.txt / ValNSIDE.txt / ConNSIDE.txt
                / DicedNSIDE.txt files. Globbing(*) or
                file names separated by " " are also
                accepted.
 -n|--name    > A string, used as a label in output
                file name and prifix of each NSIDEid,
                default is name of the first sample +
                Con.
 -o|--output  > Output folder.
 -h|--help    > Print this usage.

------- ------- ------- ------- ------- ------- -------   

USAGE
my @in_files      = '';
my $name          = '';
my $out_folder    = dirname './';
my $help          = 0;

die $usage
  unless GetOptions(
  	"i|input=s{2,}"  => \@in_files,
    "n|name:s"       => \$name,
    "o|output:s"     => \$out_folder,
    "h|help"         => \$help,
  );

##########################################################################################
# Check the parameter infomation 
##########################################################################################  
die $usage if $help;

# Check the input file
my @all_in_files;
foreach my $glob_str (@in_files) {
  push @all_in_files, grep {-s $_} glob $glob_str;
}

die "No input files\n" unless @all_in_files;
die "Multiple (>1) input files are expected\n" unless (@all_in_files >= 2);

########################################################################################## 
# Output files
##########################################################################################
# Make output directary
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $out_conNSIDE_file;

if ($name ne ''){
  $out_conNSIDE_file = $out_folder . "/ConNSIDE_" . $name .".txt";
}else{
  $out_conNSIDE_file = $out_folder . "/ConNSIDE.txt";

  my $basename =  basename($all_in_files[0], qw(.txt));
  if ($basename =~ /NSIDE_(.*)/){$name = $1;}
}

##########################################################################################
# Read Data
##########################################################################################
my %BeadData;

foreach my $in_file (sort @all_in_files){
	open( INPUT, "<", $in_file ) or die "$!";

	while(<INPUT>){
		next if (/NSIDEid/);

		chomp;
		my @F = split /\t/;

		my @Beads = split(/,/, $F[3]);
		my $num = @Beads;

		for (my $a = 0; $a < $num; $a += 1){
			my $BeadID = shift(@Beads);
			if ($BeadData{$BeadID}{LinkedBeads}){
				my @LinkedBeads = @{$BeadData{$BeadID}{LinkedBeads}};
				push @LinkedBeads, @Beads;
				$BeadData{$BeadID}{LinkedBeads} = [@LinkedBeads];
			}else{
				$BeadData{$BeadID}{LinkedBeads} = [@Beads];				
			}
			push @Beads, $BeadID;
		}
	}
	close INPUT;
}

##########################################################################################
# Join 
##########################################################################################
open( NSIDE, ">", $out_conNSIDE_file ) or die "$!";
print NSIDE "NSIDEid\tSize\tType\tBeads\n";

my $NSIDE_index = 0;

foreach my $BeadID (sort keys %BeadData){

	next if ($BeadData{$BeadID}{Visited});

	my @NSIDE = ($BeadID);

	my $visited = 0;
	my @LinkedBeads = @{$BeadData{$BeadID}{LinkedBeads}};
	$BeadData{$BeadID}{Visited} = 1;

	until ($visited == @LinkedBeads){
		$visited = 0;
		my @LinkedBeads_tmp = ();

		foreach my $LinkedBeadID (@LinkedBeads){
			if ($BeadData{$LinkedBeadID}{Visited}){
				$visited += 1;
			}else{
				push @NSIDE, $LinkedBeadID;
				push @LinkedBeads_tmp, @{$BeadData{$LinkedBeadID}{LinkedBeads}};
			}
			$BeadData{$LinkedBeadID}{Visited} = 1;
		}
		@LinkedBeads = unique(@LinkedBeads_tmp);
	}

	my $Size = @NSIDE;

	#NSIDE_type: Cis or Trans
	#NSIDE_mean_value
	my @chrs;
	my $NSIDE_type;

	my @values;

	foreach my $BeadID (@NSIDE){
		my @site = split(/:|-/, $BeadID);
		my $chr = $site[0];
		push @chrs, $chr;
	}

	if (unique(@chrs) == 1){
		$NSIDE_type = "Cis";
	}else{
		$NSIDE_type = "Trans";
	}

	#Ignore NSIDE with only one bead
	if ($Size == 1){
		next;
	}else{
		$NSIDE_index += 1;
		my $NSIDEid = $name . "_ConNISDE" . sprintf("%07d", $NSIDE_index);
		my $Beads = join ",", sort (@NSIDE);
		print NSIDE "$NSIDEid\t$Size\t$NSIDE_type\t$Beads\n";	
	}
}

close NSIDE;
