#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use PDL;
use PDL::Basic;
use PDL::Bad;
use PDL::Core;
use PDL::GSL::RNG;
use PDL::Primitive;
use Array::Utils;
use Statistics::Descriptive;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Valid.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to identify valid (repeatable)
NSIDE from replicates of 3D genome structure by local
sensitive hashing.

Reference:
1) https://www.pinecone.io/learn/locality-sensitive-
   hashing/
2) Leskovec J et al., Cambridge university press, 2020.

------- ------- ------- ------- ------- ------- -------

Iutput files (details in "Options"):

1. NSIDE.txt files of multiple replicates: identified
   NSIDE, including 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Output file:

1. ValNSIDE.txt: valid (repeatable) NSIDE of multiple 
   replicates, including at least 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.
   > If use "--detail": 
   5) Average overlap coefficient of each candidate
      NSIDE pair.
   6) Average jaccard index of each candidate NSIDE
      pair.
   7) List of NSIDEid of candidate NSIDEs.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Valid.pl -i NSIDE*.txt [Options]

NSIDE_Valid.pl -i NSIDE*.txt -nbit 100 -b 25
               -sim 3 -n name 
               -co 0.8 -cj 0.6 --detail
               -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:

 -i|--input      > NSIDE.txt files. Globbing(*) or file
                   names separated by " " are also
                   accepted.
 -nbit|--nbit    > Length of NSIDE signature (used for
                   local sensitive hashing), must be an
                   integer, default is 100. The bigger
                   nbit results in more accurate
                   simularity between NSIDEs, but
                   consumes more running time.
 -b|--band       > Number of bands (used for local
                   sensitive hashing), must be an
                   integer, default is 25. The more
                   bands results in lower false
                   negative rate in identifing
                   candidate NSIDE pairs, but consumes
                   more running time.
                 > !Notice:
                 > nbit must be evenly divisible by b!
 -sim|--simular  > Min number of structures containing
                   candidate valid NSIDE, default is
                   the number of input files.
 -n|--name       > A string, used as a label in output
                   file name and prifix of each
                   NSIDEid, default is name of the
                   first sample + Val.
 -co|--cutoffo   > Cutoff (>=) of average overlap
                   coefficient of each candidate NSIDE
                   pair, default is 0
 -cj|--cutoffj   > Cutoff (>=) of average jaccard index
                   of each candidate NSIDE pair,
                   default is 0
 --detail        > Output average overlap coefficient,
                   average jaccard index and list of
                   NSIDEid of candidate NSIDEs, default
                   is FALSE.
 -o|--output     > Output folder.
 -h|--help       > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my @in_files      = '';
my $nbit          = 100;
my $b             = 25;
my $sim           = '';
my $name          = '';
my $cutoffo       = 0;
my $cutoffj       = 0;
my $detail        = 0;
my $out_folder    = dirname './';
my $help          = 0;

die $usage
  unless GetOptions(
    "i|input=s{2,}"  => \@in_files,
    "nbit:i"    => \$nbit,
    "b|band:i"       => \$b,
    "sim|simular:i"  => \$sim,
    "n|name:s"       => \$name,
    "co|cutoffo:f"   => \$cutoffo,
    "cj|cutoffj:f"   => \$cutoffj,
    "detail"          => \$detail,
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

# Check other parameters
if ($sim){
	die "-sim is expected to be smaller than number of input files\n" unless (@all_in_files >= $sim);
}else{
	$sim = @all_in_files;
}

die "nbit must be evenly divisible by b!\n" unless ($nbit % $b == 0);

##########################################################################################
# Output Files
##########################################################################################
# Make output directary 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $out_valNSIDE_file = '';

if ($name ne ''){
  $out_valNSIDE_file = $out_folder . "/ValNSIDE_" . $name .".txt";
}else{
  $out_valNSIDE_file = $out_folder . "/ValNSIDE.txt";

  my $basename =  basename($all_in_files[0], qw(.txt));
  if ($basename =~ /NSIDE_(.*)/){$name = $1;}
}

##########################################################################################
# Read NSIDE data
##########################################################################################
my $datestring = localtime(); 
warn("$datestring: Reading NSIDE data...\n");

my %Data;
my @AllUniqBeads;

foreach my $in_file (sort @all_in_files){
	open( INPUT, "<", $in_file ) or die "$!";

	my $cell = basename($in_file, qw(.txt));
	my $cell_name = '';
	if ($cell =~ /NSIDE_(.*)/){$cell_name = $1;}

	while(<INPUT>){
	   next if (/NSIDEid/);

	   chomp;
	   my @F = split /\t/;

	   my $NSIDEid = $F[0];
	   my @Beads = split(/,/, $F[3]);
	   $Data{$cell_name}{$NSIDEid}{Beads} = [@Beads];

	   push @AllUniqBeads, @Beads;
	}
	close INPUT;
}

##########################################################################################
# Use local sensitive hashing for index of candidates
##########################################################################################
@AllUniqBeads = Array::Utils::unique(@AllUniqBeads);
my %Vocab_beads = bulid_vocab(\@AllUniqBeads);
my $arr = minhash_arr(\%Vocab_beads, $nbit);

my %buckets;
my %Vocab_cell;
my %Vocab_NSIDE;

my $cell_counter = 0;

foreach my $cell (sort keys %Data){
	$datestring = localtime(); 
	warn("$datestring: Processing $cell...\n");

	$Vocab_cell{$cell_counter} = $cell;

	my %Data_cell = %{$Data{$cell}};
	
	my @Sigs = ();
	my $NSIDE_counter = 0;
	foreach my $NSIDEid (sort keys %Data_cell){
		$Vocab_NSIDE{$cell_counter}{$NSIDE_counter} = $NSIDEid; 

		my @Sig = get_signature_list(\%Data_cell, $NSIDEid, \%Vocab_beads, $arr);
		push @Sigs, [@Sig];

		$NSIDE_counter += 1;
	}

	%buckets = add_hash($b, \@Sigs, \%buckets, $cell_counter);
	$cell_counter += 1;
}

##########################################################################################
# Verify candidates: mean overlap coefficient between candidate NSIDEs >= 0.8
##########################################################################################
$datestring = localtime(); 
warn("$datestring: Verifying candidates...\n");

my @candidates = check_candidates(\%buckets, $sim);

my %ValidData;

foreach my $candidates (@candidates){
	my @candidate = split(/,/, $candidates);

	my @overlap;
	my @jaccard;
	my @Beads_common;
	my @NSIDEs;
	for (my $a = 0; $a < @candidate - 1; $a += 1){
		my $candidate_a = $candidate[$a];
		my ($cell_idx_a, $NSIDE_idx_a) = split(/\_/,$candidate_a);
		my $cell_a = $Vocab_cell{$cell_idx_a};
		my $NSIDEid_a = $Vocab_NSIDE{$cell_idx_a}{$NSIDE_idx_a};
		my @Beads_a = @{$Data{$cell_a}{$NSIDEid_a}{Beads}};	

		for (my $b = $a + 1; $b < @candidate; $b += 1){
			my $candidate_b = $candidate[$b];
			my ($cell_idx_b, $NSIDE_idx_b) = split(/\_/,$candidate_b);
			my $cell_b = $Vocab_cell{$cell_idx_b};
			my $NSIDEid_b = $Vocab_NSIDE{$cell_idx_b}{$NSIDE_idx_b};
			my @Beads_b = @{$Data{$cell_b}{$NSIDEid_b}{Beads}};	

			my $overlap = overlap(\@Beads_a, \@Beads_a);
			my $jaccard = jaccard(\@Beads_a, \@Beads_b);
			push @overlap, $overlap;
			push @jaccard, $jaccard;

			if (@Beads_common){
				my @Beads_common_tmp = Array::Utils::unique(Array::Utils::intersect(@Beads_common, @Beads_a));
				@Beads_common = Array::Utils::unique(Array::Utils::intersect(@Beads_common_tmp, @Beads_b));
			}else{
				@Beads_common = Array::Utils::unique(Array::Utils::intersect(@Beads_a,@Beads_b)); 
			}

			push @NSIDEs, $NSIDEid_a;
			push @NSIDEs, $NSIDEid_b;			
		}
	}

	# Calculate overlap mean
	my $stat_overlap = Statistics::Descriptive::Full->new();
	$stat_overlap ->add_data(@overlap);
	my $overlap_mean = $stat_overlap->mean();

	# Calculate jaccard mean
	my $stat_jaccard = Statistics::Descriptive::Full->new();
	$stat_jaccard ->add_data(@jaccard);
	my $jaccard_mean = $stat_jaccard->mean();

	# Store results
	if ($overlap_mean >= $cutoffo && $jaccard_mean >= $cutoffj){
		my $Size = @Beads_common;

		#Ignore NSIDE with only one bead
		if ($Size > 1){
			my $Beads_common = join",", (sort @Beads_common);
			$ValidData{$Beads_common}{Overlap} = $overlap_mean;
			$ValidData{$Beads_common}{Jaccard} = $jaccard_mean;

			@NSIDEs = Array::Utils::unique(@NSIDEs);
			my $NSIDEs = join ",", (sort @NSIDEs);
			$ValidData{$Beads_common}{NSIDEs} = $NSIDEs;
		}else{
			next;
		}
	}
}

##########################################################################################
# Sort and output
##########################################################################################
$datestring = localtime(); 
warn("$datestring: Sorting and output valid NSIDEs...\n");

open( VALNSIDE, ">", $out_valNSIDE_file ) or die "$!";

if ($detail){
	print VALNSIDE "NSIDEid\tSize\tType\tBeads\tAvg_Overlap_Coefficient\tAvg_Jaccard_Index\tCandidate_NSIDEs\n";
}else{
	print VALNSIDE "NSIDEid\tSize\tType\tBeads\n";
}


my $NSIDE_index = 0;

foreach my $Beads (sort keys %ValidData){
	my @Beads = split(/,/,$Beads);
	my $Size = @Beads;

	#NSIDE_type: Cis or Trans
	#NSIDE_mean_value
	my @chrs;
	my $NSIDE_type;

	foreach my $BeadID (@Beads){
		my @site = split(/:|-/, $BeadID);
		my $chr = $site[0];
		push @chrs, $chr;
	}

	if (Array::Utils::unique(@chrs) == 1){
		$NSIDE_type = "Cis";
	}else{
		$NSIDE_type = "Trans";
	}

	$NSIDE_index += 1;
	my $NSIDE_id = $name . "_ValNSIDE" . sprintf("%07d", $NSIDE_index);

	if ($detail){
		my $overlap_mean = $ValidData{$Beads}{Overlap};
		my $jaccard_mean = $ValidData{$Beads}{Jaccard};
		my $NSIDEs = $ValidData{$Beads}{NSIDEs};
		print VALNSIDE "$NSIDE_id\t$Size\t$NSIDE_type\t$Beads\t$overlap_mean\t$jaccard_mean\t$NSIDEs\n";
	}else{
		print VALNSIDE "$NSIDE_id\t$Size\t$NSIDE_type\t$Beads\n";
	}	
}

close VALNSIDE;

##########################################################################################
# Subroutines
##########################################################################################
# Build Vocab: Generate hash of ID-index in @arr
# Parameters: \@arr
sub bulid_vocab {
	my @AllUniqBeads = @{$_[0]};
	my %Vocab;

	my $index = 0;
	foreach my $BeadID (@AllUniqBeads){
		$Vocab{$BeadID} = $index;
		$index += 1;
	}
	return %Vocab;
}

# One hot encoding: transform beads in NSIDE to a sparse vector in pdl
# Parameters: \@Beads, \%Vocab
sub one_hot {
	my @Beads = @{$_[0]};
	my %Vocab = %{$_[1]};

	my $len = %Vocab;
	my $SparseVector = zeros($len, 1);

	foreach my $BeadID (@Beads){
		my $index = $Vocab{$BeadID};
		$SparseVector -> set($index, 0, 1);
	}
	return $SparseVector;
}

# Create a shuffled hash of line numbers
# Parameters: \%Vocab, $nbit
sub minhash_arr {
	my %Vocab = %{$_[0]};
	my $nbit = $_[1];

	my $len = %Vocab;
	my $arr = zeroes($len, $nbit);

	for (my $a = 0; $a < $nbit; $a += 1){
		my $line = sequence($len);
		my $rng = PDL::GSL::RNG->new('taus');
		my $seed = $a + time();
		$rng -> set_seed($seed);
		$rng->ran_shuffle($line);
		my $row = $arr -> dice_axis(1, $a) .= $line;
	}
	return $arr;
}

# Minhash to create signature (dense vector) of the input sparse vector 
# Parameters: $arr, $SparseVector
sub get_signature {
	my $arr = $_[0];
	my $SparseVector = $_[1];

	my $Signature = minimum(index1d($arr, PDL::Primitive::which($SparseVector)));

	return $Signature;
}

# Get @Signature of specific NSIDE
# Parameters: \%Data_cell, $NSIDEid, \%Vocab, $arr
sub get_signature_list {
	my %Data_cell = %{$_[0]};
	my $NSIDEid = $_[1];
	my %Vocab = %{$_[2]};
	my $arr = $_[3];

	my @Beads = @{$Data_cell{$NSIDEid}{Beads}};
	my $SparseVector = one_hot(\@Beads, \%Vocab);
	my $Signature = get_signature($arr, $SparseVector);
	my @Signature = list($Signature);

	return @Signature;
}

# Split signature for local sensitive hashing
# Parameters: \@Sig, $b
sub make_subvecs{
	my @Sig = @{$_[0]};
	my $b = $_[1];

	my $len = @Sig;
	die "nbit must be evenly divisible by b!\n" unless ($len % $b == 0);

	my $r = int($len/$b);
	my @subvecs;

	for (my $i = 0; $i < $len; $i = $i + $r){
		my $start = $i;
		my $end = $i + $r - 1;
		my @subvec = @Sig[$start..$end];
		push @subvecs, [@subvec];
	}

	return @subvecs;
}

# Distribute subvecs of signatures to hashes
# Parameters: $b, \@Sigs, \%buckets, $cell_counter;
sub add_hash{
	my $b = $_[0];
	my @Sigs = @{$_[1]};
	my %buckets = %{$_[2]};
	my $cell_counter = $_[3];

	my $counter = 0;

	foreach my $Sig_ref (@Sigs){
		my @Sig = @{$Sig_ref};
		my @subvecs = make_subvecs(\@Sig, $b);
		my $i = 0;
		foreach my $subvec_ref (@subvecs){
			my @subvec = @{$subvec_ref};
			my $subvec = join ",", @subvec;
			if ($buckets{$i}{$subvec}){
				my @counters = @{$buckets{$i}{$subvec}};
				push @counters, "$cell_counter\_$counter";
				$buckets{$i}{$subvec} = [@counters];
			}else{
				$buckets{$i}{$subvec} = [("$cell_counter\_$counter")];
			}
			$i += 1;
		}

		$counter += 1;
	}
	return %buckets;
}

# Find candidates from built hashes
# Parameters: \%buckets, $sim
sub check_candidates{
	my %buckets = %{$_[0]};
	my $sim = $_[1];

	my @candidates;

	foreach my $i (keys %buckets){
		my @subvecs_hash_refs = $buckets{$i};
		foreach my $subvec_hash_ref (@subvecs_hash_refs){
			my %subvec = %{$subvec_hash_ref};
			foreach my $subvec (keys %subvec){
				my @counters = @{$buckets{$i}{$subvec}};
				if (@counters >= $sim){
					my $counters = join",", @counters;
					push @candidates, $counters;
				}
			}
		}
	}

	@candidates = Array::Utils::unique(@candidates);
	return @candidates;
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
