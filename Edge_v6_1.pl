#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use POSIX;
use IPC::Cmd qw[can_run run];
use List::Util qw[max min];
use threads;
use Thread::Queue;
use File::Path qw[rmtree];
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Edge.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to identify edged beads of each 
bead.

------- ------- ------- ------- ------- ------- -------

Input file:

1. RawData.txt: from NSIDE_Pre.pl, refomated and 
   annotated data of each bead (genomic bin), including 
   5 columns: BeadID, x, y, z, Value (recommand -d for 
   NSIDE_Edge.pl in the first line).

------- ------- ------- ------- ------- ------- -------

Output file:

1. EdgeBeads.txt: edged beads of each beads, including 
   2 columns: BeadID, list of BeadID of all edged beads.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Edge.pl -i RawData_xx.txt [Options]

NSIDE_Edge.pl -i RawData_xx.txt -d 2.0 -s 100 -p 16
              -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:

 -i|--input    > RawData.txt file from NSIDE_Pre.pl.
 -d|--distance > Upper limit of distance between two
                 edged beads, default 3 bead radii 
                 (inferred by NISDE_Pre.pl).
 -s|--slice    > Number of tmp file, only influence run
                 time, for 20kb structure of the whole 
                 genome, recommand -s 100 (default).
 -p|--process  > Number of threads used for calculation,
                 default 16.
 -o|--output   > Output folder.
 -h|--help     > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my $in_file           = '';
my $distance          = '';
my $slice_file_number = 100;
my $process           = 16;
my $out_folder        = dirname './';
my $help              = 0;

die $usage
  unless GetOptions(
    "i|input=s"             => \$in_file,
    "d|distance:f"          => \$distance,
    "s|slice_file_number:i" => \$slice_file_number,
    "p|process:i"           => \$process,
    "o|output:s"            => \$out_folder,
    "h|help"                => \$help,
  );

##########################################################################################
# Check the parameter infomation
##########################################################################################  
die $usage if $help;

# Check the running environment
# can_run('NSIDE_Edge_sub.pl') or die 'NSIDE_Edge_sub.pl is not executable!';

#check the input file
die "No RawData file!\n" unless -s $in_file;

##########################################################################################
# Output Files
##########################################################################################
# Make output directary 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $basename =  basename($in_file, qw(.txt));
my $name = '';
if ($basename =~ /RawData_(.*?)_(.*)/){$name = $2;}

my $out_file = '';
my $tmp_folder = '';

if ($name ne ''){
  $out_file = $out_folder . "/EdgedBeads_" . $name . ".txt";
  $tmp_folder = $out_folder . "/tmpFiles_" . $name; 
}else{
  $out_file = $out_folder . "/EdgedBeads.txt";
  $tmp_folder = $out_folder . "/tmpFiles"; 
}

# Make tmp directary 
mkdir $tmp_folder unless -d $tmp_folder;

##########################################################################################
# Divide all beads into cubes
##########################################################################################
my $distance_max;
my $split;
my $sidelength;
my $min;
my %CubeData;

open( DATA, "<", $in_file ) or die "$!";
while(<DATA>){
    next if (/BeadID/);


    if (/\A###/){
      if ($_ =~ /-d <= (.*), s = (.*), range of xyz = (.*), min of xyz = (.*), Xmean/){ $distance_max = $1, $split = $2, $sidelength = $3 / $split, $min = $4};
    }else{
      chomp;

      my @F = split /\s+/;
      my $BeadID = $F[0];
      my $x = $F[1];
      my $y = $F[2];
      my $z = $F[3];

      #divide beads into cube
      my $a = '';
      my $X = '';
      my $Y = '';
      my $Z = '';

      for ( $a = 0 ; $a < $split ; $a += 1 ){
          if($min + $a * $sidelength < $x && $x <= $min + ($a + 1) * $sidelength){$X = $a + 1;}
          if($min + $a * $sidelength < $y && $y <= $min + ($a + 1) * $sidelength){$Y = $a + 1;}
          if($min + $a * $sidelength < $z && $z <= $min + ($a + 1) * $sidelength){$Z = $a + 1;}
      }

      my $CubeID;
      if ($X && $Y && $Z){

        $CubeID = "$X-$Y-$Z";

        if ($CubeData{$CubeID}){
            push @{$CubeData{$CubeID}{BeadsArray}}, "$BeadID($x,$y,$z)";
        }else{
            $CubeData{$CubeID}{BeadsArray} = [("$BeadID($x,$y,$z)")];
        }
      }else{
        warn "Failed to divide Bead $BeadID ($x, $y, $z) into cube!\n";
      }
    }
}
close DATA;

##########################################################################################
# Generate split cube file
##########################################################################################
my @CubeIDs = sort keys %CubeData;
my $total_cube_number = @CubeIDs;
my $slice_cube_number = ceil($total_cube_number /$slice_file_number);

my @edgeBeads_cmds;

for (my $a = 1; $a <= $slice_file_number; $a += 1){
  my $tmp_slice_file = $tmp_folder . "/Slice_of_cube_$a";
  open( SLICE, ">", $tmp_slice_file) or die "$!";

  for (my $b = ($a - 1) * $slice_cube_number; $b <= ($a * $slice_cube_number) - 1; $b += 1){
    if ($CubeIDs[$b]){
      my $CubeID = $CubeIDs[$b];
      my @BeadsinCube = @{$CubeData{$CubeID}{BeadsArray}};

      my ($x, $y, $z) = split(/-/, $CubeID);
      my @EdgeCubes = ();
      my @BeadsinEdgeCubes = ();
      if ($x && $y && $z){
        for (my $x_e = $x - 1; $x_e <= $x + 1; $x_e += 1){
          for (my $y_e = $y - 1; $y_e <= $y + 1; $y_e += 1){
            for (my $z_e = $z - 1; $z_e <= $z + 1; $z_e += 1){
              if (($x_e >= 1) && ($x_e <= $split) && ($y_e >= 1) && ($y_e <= $split) && ($z_e >= 1) && ($z_e <= $split) && ("$x_e-$y_e-$z_e" ne "$x-$y-$z")){
                push @EdgeCubes, "$x_e-$y_e-$z_e";
              }else{
                next;
              }
            }
          }
        }
      }else{
          warn "Invalid CubeID: $CubeID!\n";
          next;
      }

      foreach my $EdgeCubeID (@EdgeCubes){
        if ($CubeData{$EdgeCubeID}{BeadsArray}){
          push @BeadsinEdgeCubes, @{$CubeData{$EdgeCubeID}{BeadsArray}};
        }
      }

      my $BeadsinCube = join ";", @BeadsinCube;
      my $BeadsinEdgeCubes = join ";", @BeadsinEdgeCubes;
      print SLICE "$CubeID\t$BeadsinCube\t$BeadsinEdgeCubes\n";
    }
  }
  close SLICE;

  $distance = $distance_max unless $distance;

  if (-s $tmp_slice_file){
    push @edgeBeads_cmds, "perl /data4/xyliu/Projects/NSIDE/PipelineTest/V5/Edge_sub_v5_2.pl -i $tmp_slice_file -d $distance -o $tmp_folder";  
  }
} 

run_parallel($process, @edgeBeads_cmds);

##########################################################################################
# Cat EdgeBeads files
##########################################################################################
my $cat_cmd = "cat $tmp_folder/EdgeBeads_*.txt > $out_file";
run_command($cat_cmd);

##########################################################################################
# Clean up tmp files
##########################################################################################
rmtree $tmp_folder,{verbose => 0};

##########################################################################################
# Subroutines
##########################################################################################
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
