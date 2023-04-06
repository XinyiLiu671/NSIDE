#!/usr/bin/perl -w
use strict;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Statistics::Descriptive;
use Array::Utils;
use Statistics::R;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Quantify.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to calculate jaccard index between
each NSIDE in each sample and its consensus NSIDE, or
calculate overlap coefficient of each NSIDE in each
sample and its coresponding diced NSIDEs.

------- ------- ------- ------- ------- ------- -------

Iutput files (details in "Options"):

1. NSIDE.txt, ValNSIDE.txt, ConNSIDE.txt or
   DicedNSIDE.txt files, including at least 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.
2. (Optional) ConNSIDE.txt or DicedNSIDE.txt file.

------- ------- ------- ------- ------- ------- -------

Output file:

1. Quantified_mode_SizeN_SampleM.txt: numeric matrix of
   jaccard index or overlap coefficient with cells in
   columns and consensus / diced NSIDE in rows.

> If do not use "-c":
2. ConNSIDE.txt ("consensus" mode) or DicedNSIDE.txt
   ("dice" mode) generated from input files, including
   at least 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

> If use "--heatmap":
3. Heatmap_mode_SizeN_SampleM_Clustermethod.pdf:
   (clustered) heatmap of jaccard index or overlap
   coefficient.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Quantify.pl -i NSIDE*.txt [Options]

NSIDE_Quantify.pl -i NSIDE*.txt -n name
                  -m consensus -c ConNSIDE.txt
                  -cs 10
                  --heatmap -cc 2
                  --cluster_sample --cluster_NSIDE
                  --clustering_method complete
                  --show_samplename --show_NSIDEid
                  -o output_folder_name

NSIDE_Quantify.pl -i NSIDE*.txt -n name
                  -m dice -c DicedNSIDE.txt
                  -cs 10
                  --heatmap -cc 2
                  --cluster_sample --cluster_NSIDE
                  --clustering_method ward.D
                  --show_samplename --show_NSIDEid
                  -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:
 -i|--input          > NSIDE.txt / ValNSIDE.txt /
                       ConNSIDE.txt / DicedNSIDE.txt
                       files. Globbing(*) or file names
                       separated by " " are also
                       accepted.
 -n|--name           > A string, used as a label in
                       output file name.
 -m|--mode           > Mode, "consensus" (default) or
                       "dice".
                     > In "consensus" mode: calculate
                       jaccard index between each NSIDE
                       in each sample and its consensus
                       NSIDE.
                     > In "dice" mode: calculate
                       overlap coefficient of each NSIDE
                       in each cell and its coresponding
                       diced NSIDEs.
 -c|--common_file    > (Optional) ConNSIDE.txt or
                       DicedNSIDE.txt file. If do not
                       input this file, this file will
                       be calculated from the input
                       files according to "-m".
                     > !Notice:
                     > This file must be generated from
                       exactly the same files as input
                       files. Slice of this file is also
                       accepted.
 -cs|cutoffsize      > Cutoff (>=) of consensus or
                       diced NSIDE size used to
                       calculate index and output.
                       Default is 0.
 --heatmap           > Output heatmap of jaccard index
                       or overlap coefficient. Default
                       is FALSE.
 -cc|cutoffcell      > Cutoff (>=) of sample (cell)
                       number that have index > 0 used
                       to draw the heatmap, must be an
                       integer. Default is 1.
 --cluster_sample    > Cluster samples (rows) in the
                       heatmap. Default is TRUE.
 --cluster_NSIDE     > Cluster NSIDEs (columns) in the
                       heatmap. Default is TRUE.
 --clustering_method > Clustering method used.
                       "complete" (default), "ward",
                       "ward.D", "ward.D2", "single",
                       "average", "mcquitty", "median"
                       or "centroid". Please refer to
                       the manual of R package "hclust".
 --show_samplename   > Show NSIDEid).
                       Default is FALSE.
 --show_NSIDEid      > Show column names (sample
                       names). Default is TRUE.
 -o|--output         > Output folder.
 -h|--help           > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my @in_files        = '';
my $name            = '';
my $mode            = 'consensus';
my $c_file          = '';
my $cutoffsize      = 0;
my $heatmap         = 0;
my $cutoffcell      = 1;
my $cluster_sample  = 1;
my $cluster_NSIDE   = 1;
my $method          = 'complete';
my $show_samplename = 1;
my $show_NSIDEid    = 0;
my $out_folder      = dirname './';
my $help            = 0;

die $usage
  unless GetOptions(
  	"i|input=s{2,}"       => \@in_files,
    "n|name:s"            => \$name,
  	"m|mode:s"            => \$mode,
  	"c|common_file:s"     => \$c_file,
    "cs|cutoffsize:i"     => \$cutoffsize,
    "heatmap"             => \$heatmap,
    "cc|cutoffcell:i"     => \$cutoffcell,
    "cluster_sample"      => \$cluster_sample,
    "cluster_NSIDE"       => \$cluster_NSIDE,
    "clustering_method:s" => \$method,
    "show_samplename"     => \$show_samplename,
    "show_NSIDEid"        => \$show_NSIDEid,
    "o|output:s"          => \$out_folder,
    "h|help"              => \$help,
  );

##########################################################################################
# Check the parameter infomation 
##########################################################################################  
die $usage if $help;

# Check the running environment
# can_run('NSIDE_Consensus.pl') or die 'NSIDE_Consensus.pl is not executable!';
# can_run('NSIDE_Dice.pl') or die 'NSIDE_Dice.pl is not executable!';
can_run('R') or die 'R is not executable!';
check_library('pheatmap') or die 'R package "pheatmap" has not been installed!';

# Check the input file
my @all_in_files;
foreach my $glob_str (@in_files) {
  push @all_in_files, grep {-s $_} glob $glob_str;
}

die "No input files\n" unless @all_in_files;
die "More than one files are expected\n" unless (@all_in_files >= 2);

# Check the parameters
# -m
die "-m is expected to be \"consensus\" or \"dice\"!\n" unless ($mode eq "consensus" | $mode eq "dice");

# method
my @All_methods = qw(ward ward.D ward.D2 single complete average mcquitty median centroid);
my $matched = check_paras(\@All_methods, $method);
die "Unknown clustering_method: $method!\n" unless ($matched == 1);

# If no "-c", generate common NSIDE file. 
unless ($c_file){
	if ($mode eq "consensus"){
		if ($name){
			my $consensus_cmd = "perl /data4/xyliu/Projects/NSIDE/PipelineTest/V6/Consensus_v6_1.pl -i @in_files -n $name -o $out_folder";
			run_command($consensus_cmd);
			$c_file = $out_folder . "/ConNSIDE_" . $name .".txt";			
		}else{
			my $consensus_cmd = "perl /data4/xyliu/Projects/NSIDE/PipelineTest/V6/Consensus_v6_1.pl -i @in_files -o $out_folder";
			run_command($consensus_cmd);
			$c_file = $out_folder . "/ConNSIDE.txt";						
		}
	}elsif ($mode eq "dice"){
		if ($name){
			my $dice_cmd = "perl /data4/xyliu/Projects/NSIDE/PipelineTest/V6/Dice_v6_10.pl -i @in_files -n $name -o $out_folder";
			run_command($dice_cmd);
			$c_file = $out_folder . "/DicedNSIDE_" . $name .".txt";			
		}else{
			my $dice_cmd = "perl /data4/xyliu/Projects/NSIDE/PipelineTest/V6/Dice_v6_10.pl -i @in_files -o $out_folder";
			run_command($dice_cmd);
			$c_file = $out_folder . "/DicedNSIDE.txt";	
		}
	}
}

########################################################################################## 
# Output files
##########################################################################################
# Make output directary
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $out_quantified_file;
my $out_heatmap_file;

if ($name ne ''){
  $out_quantified_file = $out_folder . "/Quantified_$mode" . "_$name" . "_Size$cutoffsize" . "_Sample$cutoffcell" . ".txt";
  $out_heatmap_file = $out_folder . "/Heatmap_$mode" . "_$name" . "_Size$cutoffsize" . "_Sample$cutoffcell" . "_Cluster$method" . ".pdf";
}else{
	$out_quantified_file = $out_folder . "/Quantified_$mode" . "_Size$cutoffsize" . "_Sample$cutoffcell" . ".txt";
  $out_heatmap_file = $out_folder . "/Heatmap_$mode" . "_Size$cutoffsize" . "_Sample$cutoffcell" . "_Cluster$method" . ".pdf";
}

##########################################################################################
# Read cNSIDE data and build vocab of cNSIDEs
##########################################################################################  
my %Data;
my %Vocab;

open( CFILE, "<", $c_file ) or die "$!";
while(<CFILE>){
	next if (/NSIDEid/);

	chomp;
	my @F = split /\t/;

	my $cNSIDEid = $F[0];
	my $Size = $F[1];

	if ($Size >= $cutoffsize){
		my @Beads = split(/,/, $F[3]);

		foreach my $BeadID (@Beads) {
			$Vocab{$BeadID} = $cNSIDEid;
		}

		$Data{$cNSIDEid} = [@Beads];
	}
}
close CFILE;

##########################################################################################
# Read input data
# "consensus" mode: calculate jaccard index between each NSIDE and its consensus NSIDE
# "dice" mode: calculate overlap coefficient of each NSIDE and coresponding diced NSIDEs
##########################################################################################  
my %Result;
my @Cells;

foreach my $in_file (sort @all_in_files){
	open( INPUT, "<", $in_file ) or die "$!";

	my $cell = basename($in_file, qw(.txt));
	my $cell_name = '';
	if ($cell =~ /NSIDE_(.*)/){$cell_name = $1;}
	push @Cells, $cell_name;

	while(<INPUT>){
	  next if (/NSIDEid/);

	  chomp;
	  my @F = split /\t/;

	  my $NSIDEid = $F[0];
	  my @Beads = split(/,/, $F[3]);

	  my $Sim;
	  if ($mode eq "consensus"){

	  	my $conNSIDEid; 	
	  	foreach my $BeadID (@Beads){
	  		if ($Vocab{$BeadID}){
	  			$conNSIDEid = $Vocab{$BeadID};
	  			last;
	  		}
	  	}

	  	if ($conNSIDEid){
			  my @conBeads = @{$Data{$conNSIDEid}};
			  
			  $Sim = jaccard(\@Beads, \@conBeads);

			  if ($Result{$conNSIDEid}{$cell_name}){
			   	my @Sims = @{$Result{$conNSIDEid}{$cell_name}};
			   	push @Sims, $Sim;
			   	$Result{$conNSIDEid}{$cell_name} = [@Sims];
			  }else{
			   	$Result{$conNSIDEid}{$cell_name} = [($Sim)];
			  }	  		
	  	}

		}elsif ($mode eq "dice"){
			my @dicedNSIDEids;
			foreach my $BeadID (@Beads){
				if ($Vocab{$BeadID}){
					push @dicedNSIDEids, $Vocab{$BeadID};
				}else{
					next;
				}
			}

			if (@dicedNSIDEids){
				@dicedNSIDEids = Array::Utils::unique(@dicedNSIDEids);

				my $dicedNSIDEids = join ",", @dicedNSIDEids;
	
				foreach my $dicedNSIDEid (@dicedNSIDEids){
					my @dicedBeads = @{$Data{$dicedNSIDEid}};
					$Sim = overlap(\@Beads, \@dicedBeads);

					if ($Result{$dicedNSIDEid}{$cell_name}){
						my @Sims = @{$Result{$dicedNSIDEid}{$cell_name}};

						print "$dicedNSIDEid\t$cell_name\t$NSIDEid\n";
						push @Sims, $Sim;
						$Result{$dicedNSIDEid}{$cell_name} = [@Sims];
					}else{
						$Result{$dicedNSIDEid}{$cell_name} = [($Sim)];
					}
				}				
			}
		}
	}
	close INPUT;
}

##########################################################################################
# Calculate mean jaccard index of each conNSIDE in each cell and output
##########################################################################################
open( QUANTI, ">", $out_quantified_file ) or die "$!";

my @colnames = ("Cells");
push @colnames, sort @Cells;
my $colnames = join "\t", @colnames;
print QUANTI "$colnames\n";

foreach my $cNSIDEid (sort keys %Result){
	my @contents = ($cNSIDEid);

	foreach my $cell_name (@Cells){
		if ($Result{$cNSIDEid}{$cell_name}){
			my @Sims =  @{$Result{$cNSIDEid}{$cell_name}};
			my $stat_sim = Statistics::Descriptive::Full->new();
			$stat_sim ->add_data(@Sims);
			my $sim_mean = $stat_sim->mean();
			push @contents, $sim_mean;			
		}else{
			push @contents, 0;
		}
	}

	my $contents = join "\t", @contents;
	print QUANTI "$contents\n";
}

##########################################################################################
# Draw heatmap
##########################################################################################
if ($heatmap){
my $R = Statistics::R->new();
my $cmds_R = <<EOF;
library(pheatmap)
data <- read.table(file = "$out_quantified_file", header = TRUE, row.names=1)
data <- t(data.frame(data))
data_row = data
data_row[data_row > 0] = 1
data <- data[which(rowSums(data_row) >= $cutoffcell),]
pheatmap(data, 
	filename = "$out_heatmap_file",
	cluster_rows = $cluster_sample,
	cluster_cols = $cluster_NSIDE,
	clustering_method = "$method",
	scale = "none", 
	show_rownames = $show_samplename,
	show_colnames = $show_NSIDEid
)
EOF

$R->run($cmds_R);
}

##########################################################################################
# Subroutines
##########################################################################################
# Check whether R package has been installed
# Parameters：$library
sub check_library {
	my $library = shift;
  my $command = "R -e \'library($library)\'";
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $command, verbose => 0 );
  return $success;
}

# Check whether the input para is one of the expected parameters
# Parameters：\@paras, $para
sub check_paras {
	my @paras = @{$_[0]};
	my $para = $_[1];

	my $matched = 0;
	foreach my $para_cur (@paras){
		if ($para_cur eq $para){
			$matched += 1;
		}
	}

	return $matched;
}

# Run commands and report errors
# Parameters：$command
sub run_command {
  my $command = shift;
  # warn "\t$command\n";
  my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
    run( command => $command, verbose => 0 );
  if ($success) {
    # warn "\tDone: $command!\n";
  }
  else {
    my @stderrs = @$stderr_buf;
    warn "\t$command:\n";
    warn "Something went wrong:\n@stderrs";
  }
}

# Calculate jaccard index (( A intersect B ) / (A union B))
# Parameters: \@arr1, \@arr2
sub jaccard {
	my @arr1 = @{$_[0]};
	my @arr2 = @{$_[1]};

	my $jaccard = Array::Utils::unique(Array::Utils::intersect(@arr1,@arr2)) / Array::Utils::unique(@arr1, @arr2);

	return $jaccard;
}

# Calculate overlap coefficent (( A intersect B ) / min(A,B))
# Parameters: \@arr1, \@arr2
sub overlap {
	my @arr1 = Array::Utils::unique(@{$_[0]});
	my @arr2 = Array::Utils::unique(@{$_[1]});

	my $overlap;
	if (@arr1 <= @arr2){
		$overlap = Array::Utils::unique(Array::Utils::intersect(@arr1,@arr2)) / @arr1;
	}else{
		$overlap = Array::Utils::unique(Array::Utils::intersect(@arr1,@arr2)) / @arr2;
	}
}
