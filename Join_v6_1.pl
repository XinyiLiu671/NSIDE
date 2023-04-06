#!/usr/bin/perl -w
use strict;
use File::Basename;
use Array::Utils qw[unique];
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Join.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to identify NSIDE by joining edged
beads with significant local Moran's I.

------- ------- ------- ------- ------- ------- -------

Iutput files:

1. SigBeads.txt: from NSIDE_LISA.pl, significant beads
   (FDR of local Moran's I < cutoff, Type = 1),
   including 10 columns:
   1)  BeadID;
   2)  Value;
   3)  Ii: local Moran's I;
   4)  zI: standardized z of Ii;
   5)  Xi: nomalized Value;
   6)  wXj: mean Xi of edged beads;
   7)  Type: 4 quadrant of Moran scatterplot, 1 for
       HHHH, 2 for HLLH, 3 for LLLL, 0 for beads with
       all edged beads value = 0, 5 for beads without
       edged beads;
   8)  pvalue of Ii;
   9)  p_bonferroni: adjusted pvalue by "Bonferroni"
       method;
   10) FDR: adjusted pvalue by "FDR" method.
2. EdgeBeads.txt: from NSIDE_Edge.pl, edged beads of
   each beads, including 2 columns: BeadID, list of
   BeadID of all edged beads.

------- ------- ------- ------- ------- ------- -------

Output file:

1. NSIDE.txt: identified NSIDE, including 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Join.pl -i SigBeads_xx.txt -e EdgedBeads_xx.txt
              [Options]

NSIDE_Join.pl -i SigBeads_xx.txt -e EdgedBeads_xx.txt
              -o output_folder_name 

------- ------- ------- ------- ------- ------- -------

Options:
 -i|--input  > SigBeads file from NSIDE_LISA.pl.
 -e|--edge   > EdgedBeads file from NSIDE_Edge.pl.
 -o|--output > Output folder.
 -h|--help   > Print this usage.

------- ------- ------- ------- ------- ------- -------   

USAGE
my $in_file       = '';
my $edge_file     = '';
my $out_folder    = dirname './';
my $help          = 0;

die $usage
  unless GetOptions(
  	"i|input=s"      => \$in_file,
    "e|edge=s"       => \$edge_file,
    "o|output:s"     => \$out_folder,
    "h|help"         => \$help,
  );

##########################################################################################
# Check the parameter infomation 
##########################################################################################  
die $usage if $help;

# Check the input file
die "No SigBeads file!\n" unless -s $in_file;
die "No EdgedBeads file!\n" unless -s $edge_file;

########################################################################################## 
# Output files
##########################################################################################
# Make output directary
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $basename =  basename($in_file, qw(.txt));
my $name = '';
if ($basename =~ /SigBeads_(.*)/){$name = $1;}

my $output_NSIDE = '';

if ($name ne ''){
  $output_NSIDE = $out_folder . "/NSIDE_" . $name . ".txt"; 
}else{
  $output_NSIDE = $out_folder . "/NSIDE.txt"; 
}

#######################################################################
# Reading SigBeads and EdgedBeads file
#######################################################################
my %SigBeadData;

open ( SIGBEAD, "<", $in_file ) or die "$!";
while(<SIGBEAD>){
	next if (/Value/);
	chomp;
	my @F = split /\s+/;
	my $SigBeadID = $F[0];
	$SigBeadData{$SigBeadID}{Value} = $F[1];
}
close SIGBEAD;

open ( EDGE, "<", $edge_file ) or die "$!";
while(<EDGE>){
	chomp;
	my @F = split /\s+/;
	my $SigBeadID = $F[0];

	if (exists $SigBeadData{$SigBeadID}){
		my @EdgedBeads = split(/,/, $F[1]);
		my @SigEdgedBeads = ();
		foreach my $EdgedBeadID (@EdgedBeads){
			if (exists $SigBeadData{$EdgedBeadID}){
				push @SigEdgedBeads, $EdgedBeadID;
			}
		}
		if (@SigEdgedBeads){
			$SigBeadData{$SigBeadID}{SigEdgedBeadsArray} = [@SigEdgedBeads];
		}else{
			delete($SigBeadData{$SigBeadID});
		}
	}

}
close EDGE;

#######################################################################
# Join SigBeads to indentify NSIDE
#######################################################################
open( NSIDE, ">", $output_NSIDE ) or die "$!";
print NSIDE "NSIDEid\tSize\tType\tBeads\n";

my $NSIDE_index = 0;

foreach my $SigBeadID (sort keys %SigBeadData){

	next if ($SigBeadData{$SigBeadID}{Visited});

	my @NSIDE = ($SigBeadID);

	my $visited = 0;
	my @SigEdgedBeads = @{$SigBeadData{$SigBeadID}{SigEdgedBeadsArray}};
	$SigBeadData{$SigBeadID}{Visited} = 1;

	until ($visited == @SigEdgedBeads){
		$visited = 0;
		my @SigEdgedBeads_tmp = ();

		foreach my $SigEdgedBeadID (@SigEdgedBeads){
			if ($SigBeadData{$SigEdgedBeadID}{Visited}){
				$visited += 1;
			}else{
				push @NSIDE, $SigEdgedBeadID;
				push @SigEdgedBeads_tmp, @{$SigBeadData{$SigEdgedBeadID}{SigEdgedBeadsArray}};
			}
			$SigBeadData{$SigEdgedBeadID}{Visited} = 1;
		}
		@SigEdgedBeads = unique(@SigEdgedBeads_tmp);
	}

	my $Size = @NSIDE;

	#NSIDE_type: Cis or Trans
	#NSIDE_mean_value
	my @chrs;
	my $NSIDE_type;

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
		my $NSIDEid = $name . "_NISDE" . sprintf("%07d", $NSIDE_index);
		my $Beads = join ",", sort (@NSIDE);
		print NSIDE "$NSIDEid\t$Size\t$NSIDE_type\t$Beads\n";		
	}
}

close NSIDE;
