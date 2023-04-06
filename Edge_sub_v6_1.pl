#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Edge_sub.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to identify edged beads in each 
slice of cube file for parallelly identifying all edged 
beads in the whole 3D genome structure.

------- ------- ------- ------- ------- ------- -------

Input file:

1. Slice_of_cube_n: beads in each cube and its edged 
cube, is tmp file from NSIDE_Edge.pl, including 3 
columns (without header): CubeID, BeadsID(x,y,z) of 
beads in this cube, BeadsID(x,y,z) of beads in edged 
cubes of this cube. 

------- ------- ------- ------- ------- ------- -------

Output file:

1. EdgeBeads.txt: edged beads of each bead, including 2
   columns: BeadID, list of BeadID of all edged beads.

------- ------- ------- ------- ------- ------- -------

Usage:
NSIDE_Edge_sub.pl -i Slice_of_cube_n [Options]

NSIDE_Edge_sub.pl -i Slice_of_cube_n -d 2.0 
                  -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:
 -i|--input    > Slice_of_cube_n file from 
                 NSIDE_Edge.pl.
 -d|--distance > Upper limit distance of edged beads,
                 default 2.0.
 -o|--output   > Output folder.
 -h|--help     > Print this usage.

------- ------- ------- ------- ------- ------- -------   

USAGE
my $in_file         = '';
my $distance        = '';
my $out_folder      = dirname './';
my $help            = 0;

die $usage
  unless GetOptions(
  	"i|input=s"       => \$in_file,
  	"d|distance:f"    => \$distance,
    "o|output:s"      => \$out_folder,
    "h|help"          => \$help,
  );

##########################################################################################
# Check the parameter infomation 
##########################################################################################  
die $usage if $help;

# Check the input file
die "No Slice of Cube file!\n" unless -s $in_file;

# Check -d
die "-d is expected!\n" unless $distance;

##########################################################################################
# Output Files
##########################################################################################
# Make output directary 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $index;
if (basename($in_file) =~ /Slice_of_cube_(.*)/){ $index = $1;}
my $out_file= $out_folder . "/EdgeBeads_" . $index . ".txt";

##########################################################################################
# For beads in each cube, scan all beads in the edged cubes
##########################################################################################
my %BeadData;

open (SLICECUBE, "<", $in_file) or die "$!";
open( OUT, ">", $out_file) or die "$!";

while(<SLICECUBE>){
    chomp;
    my @F = split /\t/;

    my $CubeID = $F[0];

    my @BeadsinCube = split(/;/, $F[1]);
    my @BeadsinEdgeCubes = split(/;/, $F[2]);
	my @AllBeads = (@BeadsinCube, @BeadsinEdgeCubes);

	foreach my $BeadData_i (sort @BeadsinCube){
		
		my @contents_i = split(/\(|,|\)/, $BeadData_i);
		my $BeadID_i = $contents_i[0];
		my @xyz_i = @contents_i[1..3];
		
		my @Edges;

		if ($BeadData{$BeadID_i}{Edge}){
			@Edges = @{$BeadData{$BeadID_i}{Edge}};
		}

		foreach my $BeadData_j (sort @AllBeads){
			my @contents_j = split(/\(|,|\)/, $BeadData_j);
			my $BeadID_j = $contents_j[0];
			my @xyz_j = @contents_j[1..3];			

			next if $BeadData{$BeadID_j}{Visited};

			if ($BeadID_j ne $BeadID_i){
				my @BeadsID = ($BeadID_i, $BeadID_j);
				my $BeadsDistance = count_beads_distance(\@xyz_i, \@xyz_j);
				if ($BeadsDistance < $distance){
					push @Edges, $BeadID_j; 
				}
			}
		}

		my $Edges = join ",", @Edges;

		print OUT "$BeadID_i\t$Edges\n";

		foreach my $BeadID (@Edges){

			next if $BeadData{$BeadID}{Visited};

			if ($BeadData{$BeadID}{Edge}){
				my @Edges = @{$BeadData{$BeadID}{Edge}};
				push @Edges, $BeadID_i;
				$BeadData{$BeadID}{Edge} = [@Edges];
				my $Edges = join ",", @Edges;
			}else{
				my @Edges;
				push @Edges, $BeadID_i;
				$BeadData{$BeadID}{Edge} = [@Edges];
				my $Edges = join ",", @Edges;
			}
		} 

		$BeadData{$BeadID_i}{Visited} = 1;
	}
}

close SLICECUBE;
close OUT;

##########################################################################################
# Clean up tmp file
##########################################################################################
unlink $in_file;

#######################################################################
# Subroutines
#######################################################################
# Calculate distance between two beads
# Parameters：\@xyz_i, \@xyz_j 
sub count_beads_distance {
	my @xyz_i = @{$_[0]};
	my @xyz_j = @{$_[1]};

	my ($x1, $y1, $z1) = @xyz_i;
	my ($x2, $y2, $z2) = @xyz_j;

	my $BeadsDistance = (($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2)**(1/2);
	return $BeadsDistance;
}
