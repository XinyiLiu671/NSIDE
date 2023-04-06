#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use IPC::Cmd qw[can_run run];
use threads;
use Thread::Queue;
use Getopt::Long;

# =cut
my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to identify NSIDE.

------- ------- ------- ------- ------- ------- -------

Iutput files (details in "Options"):

1. 3D genome structure data file(s);
2. Bed file for calculating interested value (promoter
   region of intereseted genes, CpG percentage, signal
   of histone modifications, e.g.).

------- ------- ------- ------- ------- ------- -------

Output file:

1. RawData.txt: refomated and annotated data of each
   bead (genomic bin), including 5 columns: BeadID, x, 
   y, z, Value, (recommand -d for NSIDE_Edge.pl in the 
   first line).
2. EdgeBeads.txt: from NSIDE_Edge.pl, edged beads of 
   each beads, including 2 columns: BeadID, list of 
   BeadID of all edged beads.
3. LISAresults.txt: details of local Moran's I of each
   bead, including 10 columns: 
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
4. SigBeads.txt: significant beads (FDR < cutoff, Type 
   = 1, columns are the same as LISAresults.txt).
5. NSIDE.txt: identified NSIDE, which can be further
   analyzed by NSIDE_Annotate.pl, NSIDE_Valid.pl and
   NSIDE_Diff.pl, including 4 columns: 
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE 
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE.pl -i *.3dg.txt -b Genes.1kb.promoters.bed 
             [Options]

NSIDE.pl -i *.3dg.txt -fn gene -fo 3dg_dip-c -r 20000 
         -b Promoters.bed -m count -f 1e-9
         -d 2.0 -s 100 -p 16
         -cf 0.05
         -o output_folder_name

NSIDE.pl -i *.3dg.txt -fn CpG -fo 3dg_dip-c -r 20000 
         -b CpG_Percent.bed -m direct
         -d 2.0 -s 100 -p 16
         -cf 0.05         
         -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:

 -i|--input      > 3D genome structure file(s).
                   Globbing(*) or file names separated
                   by " " are also accepted.
 -fn|--facname   > Name of factor, gene (default), 
                   H3K4me3, e.g., used as a label in
                   output file names and prifix of
                   each NSIDEid.
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
 -d|--distance   > Upper limit of distance between two
                   edged beads, default 3 bead radii 
                   (inferred by NISDE_Pre.pl).
 -s|--slice      > Number of tmp file, only influence
                   run time, for 20kb structure of the
                   whole genome, recommand -s 100 
                   (default).
 -pe|--process_e > Number of threads used for edge
                   calculation, default 16.
 -cf|--cutofffdr > Cutoff of FDR, default 0.05.
 -p|--process    > Number of threads used for
                   calculation, default 8.
 -o|--output     > Output folder.
 -h|--help       > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my @in_files          = '';
my $factor_name       = 'gene';
my $format_in         = '3dg_hickits';
my $res               = 20000;
my $bed_file          = '';
my $mode              = 'count';
my $fraction          = '1e-9';
my $distance          = '';
my $slice_file_number = 100;
my $process_e         = 16;
my $cutofffdr         = 0.05;
my $process           = 8;
my $out_folder        = dirname './';
my $help              = 0;

die $usage
  unless GetOptions(
    "i|input=s{1,}"         => \@in_files,
    "fn|facname:s"          => \$factor_name,
    "fo|format_in=s"        => \$format_in,
    "r|resolution:i"        => \$res,
    "b|bed=s"               => \$bed_file,
    "m|mode:s"              => \$mode,
    "f|fraction:f"          => \$fraction,
    "d|distance:f"          => \$distance,
    "s|slice_file_number:i" => \$slice_file_number,
    "pe|process_e:i"        => \$process_e,
    "cf|cutofffdr:f"        => \$cutofffdr,
    "p|process:i"           => \$process,
    "o|output:s"            => \$out_folder,
    "h|help"                => \$help,
  );

##########################################################################################
# Check the parameter infomation
##########################################################################################  
die $usage if $help;

# Check the running environment
# can_run('NSIDE_sub.pl') or die 'NSIDE_sub.pl is not executable!';
can_run('R') or die 'R is not executable!';

# Check the input file
my @all_in_files;
foreach my $glob_str (@in_files) {
  push @all_in_files, grep {-s $_} glob $glob_str;
}
die "No input files\n" unless @all_in_files;
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

##########################################################################################
# Get names and commands
##########################################################################################
my @NSIDE_cmds;

my @suffx = qw/.impute3.round4.3dg.txt .3dg.txt .txt/;
# my @suffx  = qw/.3dg.txt/;

foreach my $in_file (sort @all_in_files){
  my $cell_name = basename($in_file, @suffx);
  my $cell_out_folder = $out_folder . "/" . $cell_name;
  my $cmd = "perl /data4/xyliu/Projects/NSIDE/PipelineTest/V6/NSIDE_sub_v5_1.pl -i $in_file -fn $factor_name -n $cell_name -fo $format_in -r $res -b $bed_file -f $fraction -m $mode -d $distance -s $slice_file_number -p $process_e -cf $cutofffdr -o $cell_out_folder";
  push @NSIDE_cmds, $cmd;
}

run_parallel($process, @NSIDE_cmds);

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

# Given commands, run them in multiple threads
# Parameters：$process_num, @commonds
sub run_parallel{
  my ($process_num, @missions) = @_;
  my $stream = Thread::Queue->new(@missions,undef);
  my $mission_num = scalar @missions;

  #assgn the task
  my @running = ();
  my @Threads;
  while (@Threads < @missions) {
      @running = threads->list(threads::running);

      if ( @running < $process_num) {
      my $command = $stream->dequeue();
          my $thread = threads->new(\&run_command,$command);
          push (@Threads, $thread);
          my $tid = $thread->tid;
      }
      @running = threads->list(threads::running);
      foreach my $thr (@Threads) {
          if ($thr->is_running()) {
                my $tid = $thr->tid;
          }
          elsif ($thr->is_joinable()) {
                my $tid = $thr->tid;
                $thr->join;
          }
      }
 
      @running = threads->list(threads::running);
  }

  #join the threads
  while (@running) {
        foreach my $thr (@Threads) {
            $thr->join if ($thr->is_joinable());
        }
        @running = threads->list(threads::running);
        sleep(3);
  }
  return 0;
}
