#!/usr/bin/perl -w
use strict;
use File::Basename;
use Array::Utils;
use List::Util qw[pairs];
use List::MoreUtils;
use PDL;
use PDL::Bad;
use PDL::Core;
use Data::Dumper;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Consensus.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to dice NSIDEs into smallest sets
of beads that always be in a single NSIDE or disappear
togther in two or more cells (in other words, finding
covarying beads).

------- ------- ------- ------- ------- ------- -------

Iutput files:

1. NSIDE.txt, ValNSIDE.txt, ConNSIDE.txt or
   DicedNSIDE.txt files, including at least 4 columns:
   including at least 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Output file:

1. DicedNSIDE.txt: consensus NSIDE, including 4
   columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_Dice.pl -i NSIDE*.txt [Options]

NSIDE_Dice.pl -i NSIDE*.txt -n name
              --keepone
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
                Diced.
 --keepone    > By defaul, this program will dice all
                cell-specific NSIDE (that show up in
                only one cell) into single bead. Please
                use "--keepone" to avoid this.
 -o|--output  > Output folder.
 -h|--help    > Print this usage.

------- ------- ------- ------- ------- ------- -------   

USAGE
my @in_files      = '';
my $name          = '';
my $keepone       = 0;
my $out_folder    = dirname './';
my $help          = 0;

die $usage
  unless GetOptions(
  	"i|input=s{2,}"  => \@in_files,
    "n|name:s"       => \$name,
    "keepone"        => \$keepone,
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

$PDL::whichND = 's';

########################################################################################## 
# Output files
##########################################################################################
# Make output directary
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $out_dicedNSIDE_file;

if ($name ne ''){
  $out_dicedNSIDE_file = $out_folder . "/DicedNSIDE_" . $name .".txt";
}else{
  $out_dicedNSIDE_file = $out_folder . "/DicedNSIDE.txt";

  my $basename =  basename($all_in_files[0], qw(.txt));
  if ($basename =~ /NSIDE_(.*)/){$name = $1;}
}

##########################################################################################
# Read input data
##########################################################################################
my @AllBeads;

my $NSIDE_counter = 0;
foreach my $in_file (sort @all_in_files){
	open( INPUT, "<", $in_file ) or die "$!";

	while(<INPUT>){
	   next if (/NSIDEid/);

	   chomp;
	   my @F = split /\t/;

	   my @Beads = split(/,/, $F[3]);

	   push @AllBeads, [@Beads];

	   $NSIDE_counter += 1;
	}
	close INPUT;
}

##########################################################################################
# Dice NSIDE
##########################################################################################
# Get all uniq beads
my @AllUniqBeads;
foreach my $ref (@AllBeads){
	push @AllUniqBeads, @{$ref};
}
@AllUniqBeads = sort(Array::Utils::unique(@AllUniqBeads));

# Get pdl array of one-hot sparse vectors
my %Vocab_beads = bulid_vocab(\@AllUniqBeads);
my $len = @AllUniqBeads; 
my $SparseVectors = zeroes($len, $NSIDE_counter);
my $counter = 0;
foreach (@AllBeads){
	my @Beads = @{$_};
	my $SparseVector = one_hot(\@Beads, \%Vocab_beads);
	my $row = $SparseVectors -> dice_axis(1, $counter) .= $SparseVector;
	$counter += 1;
}

# Dice
## For common and specific
my $FinalResults_idxes_ref;
my $Col_idxes_known_ref;
if ($keepone){
	($FinalResults_idxes_ref, $Col_idxes_known_ref) = common_specific_keepone($SparseVectors);
}else{
	($FinalResults_idxes_ref, $Col_idxes_known_ref) = common_specific($SparseVectors);
}
my @FinalResults_idxes = @{$FinalResults_idxes_ref};
my @Col_idxes_known = @{$Col_idxes_known_ref};
my @Col_idxes_all = (0..$#AllUniqBeads);

## Non-common and non-specific
my @Col_idxes_unknown = Array::Utils::array_minus(@Col_idxes_all, @Col_idxes_known);
@Col_idxes_unknown = grep defined, @Col_idxes_unknown;

if (@Col_idxes_unknown){
	# Cut out known columns
	my $Col_idxes_unknown = pdl(@Col_idxes_unknown);
	my $matrix = $SparseVectors -> dice($Col_idxes_unknown);

	# Get all col sum
	my $Col_sum = sumover(transpose($matrix));
	my @Col_sum_uniq = Array::Utils::unique(list($Col_sum));
	@Col_sum_uniq = sort({$b<=>$a} @Col_sum_uniq);

	foreach my $target_Col_sum (@Col_sum_uniq){
		# Refresh unknown, matrix and col_sum
		$Col_idxes_unknown = pdl(@Col_idxes_unknown);
		$matrix = $SparseVectors -> dice($Col_idxes_unknown);
		$Col_sum = sumover(transpose($matrix));

		my $Col_idxes_diced = which($Col_sum == $target_Col_sum);
		my @Col_idxes_diced = list($Col_idxes_diced);

		# Get all row sum
		my $Row_sum = sumover($matrix -> dice($Col_idxes_diced));
		my @Row_sum_uniq = Array::Utils::unique(list($Row_sum));
		my $Col_idxes_diced_len = @Col_idxes_diced;
		my @arr_pretreat = ($Col_idxes_diced_len-1, 0);
		@Row_sum_uniq = Array::Utils::array_minus(@Row_sum_uniq, @arr_pretreat);
		@Row_sum_uniq = sort({$a<=>$b} @Row_sum_uniq);

		# Decare variables
		my $target_Row_sum;
		my $Row_idxes;
		my $matrix_diced;
		my @Col_idxes_known_thisCS;

		##> $target_Row_sum == @Col_idxes_diced - 1
		$target_Row_sum = $Col_idxes_diced_len - 1;
		$Row_idxes = which($Row_sum == $target_Row_sum);
		if (any $Row_idxes >= 0){
			$matrix_diced = $matrix -> dice($Col_idxes_diced, $Row_idxes);

			my @position_zero = list(transpose(whichND($matrix_diced==0)->dice(0)));
			my @Col_idxes_known_cur = @Col_idxes_unknown[@Col_idxes_diced[@position_zero]];

			push @Col_idxes_known, @Col_idxes_known_cur;
			push @Col_idxes_known_thisCS, @Col_idxes_known_cur;
			push @FinalResults_idxes, @Col_idxes_known_cur;
		}
		next if (@Col_idxes_known_thisCS == @Col_idxes_diced);

		##> Other $target_Row_sum: foreach rows from min to max
		$target_Row_sum = shift @Row_sum_uniq;
		if ($target_Row_sum){
			$Row_idxes = which($Row_sum == $target_Row_sum);
			$matrix_diced = $matrix -> dice($Col_idxes_diced, $Row_idxes);
			my @nested_idxes = get_nested_idxes($matrix_diced, $target_Row_sum, \@Col_idxes_unknown, \@Col_idxes_diced);

			foreach my $target_Row_sum_cur (@Row_sum_uniq){
				$Row_idxes = which($Row_sum == $target_Row_sum_cur);
				$matrix_diced = $matrix -> dice($Col_idxes_diced, $Row_idxes);

				my @nested_idxes_cur = get_nested_idxes($matrix_diced, $target_Row_sum_cur, \@Col_idxes_unknown, \@Col_idxes_diced);
				@nested_idxes = same_diff_arr_nested(\@nested_idxes, \@nested_idxes_cur);
			}

			foreach (@nested_idxes){
				my @Col_idxes_known_cur = @{$_};
				push @Col_idxes_known, @Col_idxes_known_cur;
				push @Col_idxes_known_thisCS, @Col_idxes_known_cur;
				push @FinalResults_idxes, (join ",", (sort({$a<=>$b} @Col_idxes_known_cur)));
			}	
			next if (@Col_idxes_known_thisCS == @Col_idxes_diced);				
		}
		
		# Refresh unknown
		@Col_idxes_unknown = Array::Utils::array_minus(@Col_idxes_all, @Col_idxes_known);
		@Col_idxes_unknown = grep defined, @Col_idxes_unknown;
		last unless (@Col_idxes_unknown);
	}
}

##########################################################################################
# Sort and output
##########################################################################################
my @FinalResults;

foreach my $Beads_idxes (@FinalResults_idxes){
	my @Beads_idxes = split(/,/, $Beads_idxes);

	my @Result;
	foreach my $Bead_idx (@Beads_idxes){
		push @Result, $AllUniqBeads[$Bead_idx];
	}
	my $Result = join ",", @Result;

	push @FinalResults, $Result;
}


open( NSIDE, ">", $out_dicedNSIDE_file ) or die "$!";
print NSIDE "NSIDEid\tSize\tType\tBeads\n";

my $NSIDE_index = 0;
foreach my $Beads (sort @FinalResults){
	my @NSIDE = split(/,/, $Beads);

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

	if (Array::Utils::unique(@chrs) == 1){
		$NSIDE_type = "Cis";
	}else{
		$NSIDE_type = "Trans";
	}

	$NSIDE_index += 1;
	my $NSIDEid = $name . "_DicedNISDE" . sprintf("%07d", $NSIDE_index);
	my $Beads = join ",", sort (@NSIDE);
	print NSIDE "$NSIDEid\t$Size\t$NSIDE_type\t$Beads\n";
}

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

# Get combos and array of col index of common and specific elements in a one-hot pdl matrix
# Parameter: $matrix
sub common_specific {
	my $matrix = $_[0];

	my $Col_sum = sumover(transpose($matrix));
	my $Row_num = (dims($matrix))[1];

	my @FinalResults_idxes;
	my @Col_idxes_known;

	## sum == total NSIDE number: all in one set
	my @Col_idxes = list(which($Col_sum == $Row_num));
	my @Result_idx;
	if (@Col_idxes){
		foreach my $Col_idx (@Col_idxes){
			push @Result_idx, $Col_idx;
		}
		push @FinalResults_idxes, (join ",", (sort({$a<=>$b} @Result_idx)));		
	}
	push @Col_idxes_known, @Col_idxes;

	## sum == 1: each in one set
	@Col_idxes = list(which($Col_sum == 1));
	if (@Col_idxes){
		foreach my $Col_idx (@Col_idxes){
			push @FinalResults_idxes, $Col_idx;
		}
	}
	push @Col_idxes_known, @Col_idxes;

	return (\@FinalResults_idxes, \@Col_idxes_known);
}

# Get combos and array of col index of common and specific elements in a one-hot pdl matrix
# Parameter: $matrix
sub common_specific_keepone {
	my $matrix = $_[0];

	my $Col_sum = sumover(transpose($matrix));
	my $Row_num = (dims($matrix))[1];

	my @FinalResults_idxes;
	my @Col_idxes_known;
	
	my @Result_idx;

	## sum == total NSIDE number: all in one set
	my @Col_idxes = list(which($Col_sum == $Row_num));
	if (@Col_idxes){
		foreach my $Col_idx (@Col_idxes){
			push @Result_idx, $Col_idx;
		}
		push @FinalResults_idxes, (join ",", (sort({$a<=>$b} @Result_idx)));		
	}
	push @Col_idxes_known, @Col_idxes;

	## sum == 1: each row in one set
	@Col_idxes = list(which($Col_sum == 1));
	if (@Col_idxes){
		my $matrix_dice = $matrix-> dice(which($Col_sum == 1));
		my @position_one = list(whichND($matrix_dice==1));

		my %Result;
		my $x;
		my $y;
		foreach my $xy (pairs @position_one){
			($x, $y) = @{$xy};

			if ($Result{$y}){
				my @xs = @{$Result{$y}};
				push @xs, $x;
				$Result{$y} = [@xs];
			}else{
				$Result{$y} = [($x)];
			}
		}

		foreach my $y_cur (keys %Result){
			@Result_idx = @{$Result{$y_cur}};
			push @FinalResults_idxes, (join ",", (sort({$a<=>$b} @Result_idx)));
		}
	
	}
	push @Col_idxes_known, @Col_idxes;

	return (\@FinalResults_idxes, \@Col_idxes_known);
}

# Dice an nested array by an flatten array
# Parameter: \@input_arr, \@para_arr
sub same_diff_arr {
	my @input_arr = @{$_[0]};
	my @para_arr = @{$_[1]};

	my @output_arr;

	foreach my $ref (@input_arr){

		if (@para_arr){
			my @input_arr_cur = @{$ref};
			my @intersect = Array::Utils::intersect(@input_arr_cur, @para_arr);

			if (@intersect){
				push @output_arr, [@intersect];

				my @minus = Array::Utils::array_minus(@input_arr_cur, @para_arr);
				if (@minus){
					push @output_arr, [@minus];
				}

				@para_arr = Array::Utils::array_minus(@para_arr, @intersect);

			}else{
				push @output_arr, [@input_arr_cur];
			}

		}else{
			push @output_arr, $ref;
		}
	}

	if (@para_arr){
		push @output_arr, [@para_arr];
	}

	return @output_arr;
}

# Dice an nested array by an nested array
# Parameter: \@arr, \@para_arr
sub same_diff_arr_nested {
	my @arr = @{$_[0]};
	my @para_arr = @{$_[1]};

	foreach my $ref (@para_arr){
		my @para_arr_cur = @{$ref};
		@arr = same_diff_arr(\@arr, \@para_arr_cur);

	}
	return @arr;
}

# Generate a nested array of all uniq combos of one in a one-hot pdl matrix
# Parameter: $matrix, $one_num, \@Col_idxes_unknown, \@Col_idxes_diced
sub get_nested_idxes {
	my $matrix_diced = $_[0];
	my $one_num = $_[1];
	my @Col_idxes_unknown = @{$_[2]};
	my @Col_idxes_diced = @{$_[3]};

	my @FinalResults_idxes_nested;

	my @position_one = list(transpose(whichND($matrix_diced==1)->dice(0)));

	if ((dims($matrix_diced))[1] == 1){
	  my @Col_idxes_cur = @Col_idxes_unknown[@Col_idxes_diced[@position_one]];
	  @FinalResults_idxes_nested = ([@Col_idxes_cur]);

	}elsif ((dims($matrix_diced))[1] >= 2){
		my $first_end = $one_num - 1;
		my @tmp_Result_idx_first =  @Col_idxes_unknown[@Col_idxes_diced[@position_one[0..$first_end]]];
		my @arr = ([@tmp_Result_idx_first]);

		for (my $i = $first_end + 1; $i < (@position_one - $first_end); $i += $one_num){
			my $end = $i + $one_num - 1;
			my @tmp_Result_idx_cur =  @Col_idxes_unknown[@Col_idxes_diced[@position_one[$i..$end]]];
			@arr = same_diff_arr(\@arr, \@tmp_Result_idx_cur);
		}

		@FinalResults_idxes_nested = @arr;
	}

	return @FinalResults_idxes_nested;
}
