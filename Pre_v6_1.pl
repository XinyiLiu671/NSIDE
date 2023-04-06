#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Statistics::Descriptive;
use List::Util qw[max min];
use POSIX;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Pre.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to reformat 3D genome structure 
data in various formats from hickits(.3dg), dip-c(.3dg.
txt), nuc_dynamics(.pdb), Chrom3D(.cmm), MERFISH, 
seqFISH or data in xyz format and annotate genomic bins 
with interested values (counts of gene number, CpG 
percentage, signal of histone modifications, e.g.) for 
NSIDE identification.

------- ------- ------- ------- ------- ------- -------

Iutput files (details in "Options"):

1. 3D genome structure data file;
2. Bed file for calculating interested value (promoter 
   region of intereseted genes, CpG percentage, signal 
   of histone modifications, e.g.).

------- ------- ------- ------- ------- ------- -------

Output file:

1. RawData.txt: refomated and annotated data of each 
   bead (genomic bin), including 5 columns: BeadID, x, 
   y, z, Value, (recommand -d for NSIDE_Edge.pl in the 
   first line).

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Pre.pl -i XXX.3dg.txt -b Genes.1kb.promoters.bed 
             [Options]

NSIDE_Pre.pl -i XXX.3dg.txt -n Cell1 -fn gene
             -fo 3dg_dip-c -r 20000 
             -b Promoters.bed -m count -f 1e-9
             -o output_folder_name

NSIDE_Pre.pl -i XXX.3dg.txt -n Cell1 -fn CpG
             -fo 3dg_dip-c -r 20000 
             -b CpG_Percent.bed -m direct
             -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:

 -i|--input      > 3D genome structure file.
 -fn|--facname   > Name of factor, gene (default),
                   H3K4me3, e.g., used as a label in
                   output file name.
 -n|--name       > Name of sample, used as a label in
                   output file name.
 -fo|--format    > Format of input file, "3dg_hickits"
                   (default), "3dg_dip-c" , "pdb",
                   "cmm" or "xyz" (with genomic region
                   chr:start-end in the 5th column, if
                   not, annotate each bead with pseudo
                   value = 1).
 -r|--resolution > Resolution (in bp) of input file,
                   default 20000 (20kb).
 -b|--bed        > Bed file for calculating interested
                   value (promoter region of
                   intereseted genes, CpG percentage,
                   signal of histone modifications,
                   e.g.).
                 > !Notice:
                 > For "-m direct", bed file must be
                   of the same reference genome version
                   and resolution as the input file
                   with values in the 4th column.
 -m|--mode       > Mode of annotation, "count"
                   (default) or "direct":
                 > "count": for each bin in the input
                   file count the number of hits in the
                   bed file, while restricting to -f;
                 > "direct": directly annotate each bin
                   according to the 4th column of bed 
                   file.
 -f|--fraction   > Minimum overlap required as a
                   fraction of resolution, defalut 1e-9
                   (i.e. 1bp).
 -o|--output     > Output folder.
 -h|--help       > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my $in_file       = '';
my $factor_name   = 'gene';
my $format_in     = '3dg_hickits';
my $res           = 20000;
my $bed_file      = '';
my $mode          = 'count';
my $fraction      = '1e-9';
my $name          = '';
my $out_folder    = dirname './';
my $help          = 0;

die $usage
  unless GetOptions(
    "i|input=s"      => \$in_file,
    "fn|facname:s"   => \$factor_name,
    "fo|format_in=s" => \$format_in,
    "r|resolution:i" => \$res,
    "b|bed=s"        => \$bed_file,
    "m|mode:s"       => \$mode,
    "f|fraction:f"   => \$fraction,     
    "n|name:s"       => \$name,
    "o|output:s"     => \$out_folder,
    "h|help"         => \$help,
  );

##########################################################################################
# Check the parameter infomation
##########################################################################################  
die $usage if $help;

# Check the running environment
can_run('bedtools') or die 'bedtools is not installed!';

# Check the input file
die "No input file!\n" unless -s $in_file;
die "No bed file!\n" unless -s $bed_file;

# Factor name: change "_" to "-"
$factor_name =~ s/_/-/;

# Check the parameters
# -fo
my @All_formats = qw(3dg_hickits 3dg_dip-c pdb cmm xyz);
my $matched = check_paras(\@All_formats, $format_in);
die "Unknown format of input file: $format_in!\n" unless ($matched == 1);

# -m
die "-m is expected to be \"count\" or \"direct\"!\n" unless ($mode eq "count" | $mode eq "direct");

##########################################################################################
# Output Files
##########################################################################################
# Make output directary 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $out_rawdata_file = '';
my $tmp_beadsbed_file = '';

if ($name ne ''){
  $out_rawdata_file = $out_folder . "/RawData_$factor_name" . "_$name.txt";
  $tmp_beadsbed_file = $out_folder . "/tmp_beads_$factor_name" . "_$name.bed";
}else{
  $out_rawdata_file = $out_folder . "/RawData_$factor_name" . ".txt";
  $tmp_beadsbed_file = $out_folder . "/tmp_beads_$factor_name" . ".bed";
}

##########################################################################################
# Pre-process file to generate tmp bed file 
##########################################################################################
# Outputing bed of beads
open( INPUT, "<", $in_file ) or die "$!";
open( TMPBEADS, ">", $tmp_beadsbed_file) or die "$!";

if ($format_in eq "3dg_dip-c"){
  while(<INPUT>){
    chomp;
    my @F = split /\s+/;
    my @chr = split(/\(|-|_|\)/, $F[0]);
    (my $chr = $chr[0] ) =~ s/chr//g;
    my $start = $F[1];
    my $end = $F[1] + $res;
    my $xyz = join "\t", @F[2..4];
    if ($chr[1]){
      print TMPBEADS "chr$chr\t$start\t$end\t$xyz\tchr$chr\_$chr[1]:$start-$end\n";   
    }else{
      print TMPBEADS "chr$chr\t$start\t$end\t$xyz\tchr$chr:$start-$end\n";       
    }
  }
}

close INPUT;
close TMPBEADS;

##########################################################################################
# Generate RawData file
##########################################################################################
#Annotate beads with numbers of interested genomic elements or other values
my $intersect_cmd;

if ($format_in ne "xyz"){
  

  if ($mode eq "count"){
    $intersect_cmd = "intersectBed -a $tmp_beadsbed_file -b $bed_file -c -f $fraction | perl -lane '\@site = split(/:|-/, \$F[6]); \$site = join \"\\t\", \@site; \$content = join \"\\t\", \@F[3..7]; print \"\$site\\t\$content\"' | sortBed -i -| awk '{print \$7\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$8}' | sed -e \"1iBeadID\\tx\\ty\\tz\\tValue\" > $out_rawdata_file";
  }elsif ($mode eq "direct"){
    $intersect_cmd = "intersectBed -a $tmp_beadsbed_file -b $bed_file -f $fraction -wo | perl -lane '\@site = split(/:|-/, \$F[6]); \$site = join \"\\t\", \@site; \$value =  sprintf(\"%.0f\",\$F[10]); \$content = join \"\\t\", (\@F[3..6], \$value); print \"\$site\\t\$content\"' | sortBed -i -| awk '{print \$7\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$8}' | sed -e \"1iBeadID\\tx\\ty\\tz\\tValue\" > $out_rawdata_file";
    # $intersect_cmd = "intersectBed -a $tmp_beadsbed_file -b $bed_file -c -f $fraction -wo > $out_rawdata_file";
  }
}

run_command($intersect_cmd); 

##########################################################################################
# Recommand -d fot NSIDE_Edge.pl:
# According to distribution of distance between adjacent beads (~2 infered bead radii)
# 95% of adjacent beads must be considered as edged beads
# For -d: d(adjacent beads) ~ 2r, d>2r && d<=3r

# Caluculate max $split for NSIDE_Edge.pl:  
# sidelength (of split cube) > d
##########################################################################################
my $chr_pre;
my $x_pre;
my $y_pre;
my $z_pre;
my $end_pre;

my @d_adj;

my @x;
my @y;
my @z;
my @value;

#Reading data
open( DATA, "<", $out_rawdata_file ) or die "$!";
while(<DATA>){
    next if (/BeadID/);
    
    chomp;

    my @F = split /\s+/;
    my @id = split (/:|-/, $F[0]);
    my $chr = $id[0];
    my $start = $id[1];
    my $end = $id[2];
    my $x = $F[1];
    my $y = $F[2];
    my $z = $F[3];
    my $value = $F[4];
    # my $value = sprintf("%.0f",$F[4]);
    my $d;

    push @x, $x;
    push @y, $y;
    push @z, $z;
    push @value, $value;

    if ($chr_pre){
      if (($chr eq $chr_pre)&&($start == $end_pre)){
        $d = (($x-$x_pre)**2 + ($y-$y_pre)**2 + ($z-$z_pre)**2)**(1/2);
        push @d_adj, $d;

        $chr_pre = $chr;
        $end_pre = $end;
        $x_pre = $x;
        $y_pre = $y;
        $z_pre = $z;   
      }else{
        $chr_pre = $chr;
        $end_pre = $end;
        $x_pre = $x;
        $y_pre = $y;
        $z_pre = $z;        
      }
    }else{
      $chr_pre = $chr;
      $end_pre = $end;
      $x_pre = $x;
      $y_pre = $y;
      $z_pre = $z;
    }
}
close DATA;

# Recommand -d and calculate max split, Xmean and Xsd for LISA
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@d_adj);
my $d_adj_p95 = $stat->percentile(95);
my $min = floor(min(min(@x), min(@y), min(@z)));
my $range  = ceil(max(max(@x), max(@y), max(@z))) - $min;

my $max_split = floor($range / $d_adj_p95);

my $min_d = sprintf("%.2f", $d_adj_p95);
my $max_d = sprintf("%.2f", 1.5*$d_adj_p95);

my $stat2 = Statistics::Descriptive::Full->new();
$stat2->add_data(@value);
my $Xmean = $stat2->mean();
my $Xsd = $stat2->standard_deviation();

my $sed_cmd = "sed -i \"1i### $min_d < -d <= $max_d, s = $max_split, range of xyz = $range, min of xyz = $min, Xmean = $Xmean, Xsd = $Xsd\" $out_rawdata_file";
run_command($sed_cmd); 

##########################################################################################
# Clean up tmp file
##########################################################################################
unlink $tmp_beadsbed_file;

##########################################################################################
# Subroutines
##########################################################################################
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
