#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use IPC::Cmd qw[can_run run];
use List::MoreUtils qw[uniq];
use Statistics::Descriptive;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_Annotate.pl V5.17, written by Xinyi Liu
------- ------- ------- ------- ------- ------- -------

This program is used to annotate NSIDE by genes, other
genomic elements, or any other values of each genomic
bin.

------- ------- ------- ------- ------- ------- -------

Iutput files (details in "Options"):

1. NSIDE.txt: identified NSIDE, including 4 columns:
   1) NSIDEid;
   2) Size of NSIDE: number of beads in the NSIDE;
   3) Type of NSIDE: "Cis" means all beads in the NSIDE
      belong to same chromatin, otherwise "Trans";
   4) List of BeadID of all beads in the NSIDE.
2. Bed file for annotation (promoter region of
   intereseted genes, signal of histone modifications,
   e.g.), with id of elements or values in the
   specified column(s).
3. (Optional) In "element" mode, annotation file
   containing interested values of elements in the
   specified column(s) (e.g. gene expression level) is
   accepted.

------- ------- ------- ------- ------- ------- -------

Output file:

1. NSIDE_annotation.txt: annotated NSIDE, including at
   least 5 columns:
   > In "element" mode:
     1-4) Same as NSIDE.txt;
     5)   Number of elements;
     6-n) If input annotation file: specified
          statistics of values (by --method) in each
          specfied column (by -c).
     n+1) If use "--detail": list of elements.
   > In "direct" mode:
     1-4) Same as NSIDE.txt;
     5-n) Specified statistics of values (by --method)
          in each specfied column (by -c) of the input 
          bed file.

------- ------- ------- ------- ------- ------- -------

Usage:
NSIDE_Annotate.pl -i NSIDE.txt -b Promoters.bed
                  [Options]

NSIDE_Annotate.pl -i NSIDE.txt -b Promoters.bed
                  --header_bed
                  -m element -f 1e-9
                  -a Annotation.txt -c 5,1,7
                  --method mean,median,range
                  --header_annotate
                  --detail
                  -o output_folder_name

NSIDE_Annotate.pl -i NSIDE.txt -b CpG_Percentage.bed
                  --header_bed
                  -m direct -c 10:7
                  --method min,max,std
                  -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:

 -i|--input        > NSIDE.txt file from NSIDE_Join.pl.
 -b|--bed          > Bed file for annotation (promoter
                     region of intereseted genes,
                     signal of histone modifications,
                     e.g.), with id of elements or
                     values in the specified column(s).
                   > !Notice:
                   > For "-m element", ids of elements
                     are expected in the 4th column.
                   > For "-m direct", bed file must be
                     of the same reference genome
                     version and resolution as the
                     input with values in the 4th (more
                     columns are allowed).
                   > Header is also allowed and in 
                     "direct" mode the column name of 
                     each annotation column will be 
                     showed in the output file (please
                     use --header_bed if the file 
                     contains header, default is FALSE).
 --header_bed      > Use this option if the bed file
                     contains header, default is FALSE.
 -m|--mode         > Mode of annotation, "element"
                     (default) or "direct":
                   > "element": for each bin in the
                     input file, annotate ids of all
                     overlapping elements in the bed
                     file, while restricting to -f;
                   > "direct": directly annotate each
                     bin according to the 4th (or more)
                     column(s) of bed file, which can
                     be specified by -c.
 -f|--fraction     > Minimum overlap required as a
                     fraction of resolution, defalut
                     1e-9 (i.e. 1bp).
 -a|--annotate     > For "element" mode, annotation
                     file containing interested values
                     of elements in the specified
                     column(s) (by -c) is accepted.
                     Columns in this file should be
                     separated by tab.
                   > Header is also allowed and the
                     column name of each annotation
                     column will be showed in the
                     output file (please use
                     --header_annotate if the file
                     contains header, default is TRUE).
                   > !Notice:
                     Element id must be in the 1st
                     column and consistent to the input
                     bed file.
 --header_annotate > Use this option if the annotation
                     file contains header, default is
                     TRUE.
 -c|--column       > 1-based column number, such as
                     "5,17,10:7" is accepted, "," is
                     for separating and ":" means a
                     range of columns, default is 2 for
                     "element" mode and 4 for "direct"
                     mode.
 --method          > "mean"(default), "median", "std",
                     "range", "min", "max". More than
                     one method separated by "," is 
                     accepted, such as "mean,std,range".
 --detail          > In "element" mode, output list of
                     elements in each NSIDE in the
                     output file, default is FALSE.
 -o|--output       > Output folder.
 -h|--help         > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my $in_file         = '';
my $bed_file        = '';
my $header_bed      = 0;
my $mode            = 'element';
my $fraction        = '1e-9';
my $annotate_file   = '';
my $header_annotate = 1;
my $column_str      = '';
my $method          = 'mean';
my $detail          = 0;
my $out_folder      = dirname './';
my $help            = 0;

die $usage
  unless GetOptions(
    "i|input=s"       => \$in_file,
    "b|bed=s"         => \$bed_file,
    "header_bed"      => \$header_bed,   
    "m|mode:s"        => \$mode,
    "f|fraction:f"    => \$fraction,
    "a|annotate:s"    => \$annotate_file,
    "header_annotate" => \$header_annotate,
    "c|columns:s"     => \$column_str,
    "method:s"        => \$method,
    "detail"          => \$detail,
    "o|output:s"      => \$out_folder,
    "h|help"          => \$help,
  );

##########################################################################################
# Check the parameter infomation
##########################################################################################  
die $usage if $help;

# Check the running environment
can_run('bedtools') or die 'bedtools is not installed!';

# Check the input file
die "No NSIDE file!\n" unless -s $in_file;
die "No bed file!\n" unless -s $bed_file;

# Check the parameters
# -m
die "-m is expected to be \"element\" or \"direct\"!\n" unless ($mode eq "element" | $mode eq "direct");

# -c
unless ($column_str){
	if ($mode eq "element"){
		$column_str = "2";
	}elsif ($mode eq "direct"){
		$column_str = "4";
	}
}
my @column_indexs = get_columns($column_str);

# --method
my @method = split(/,/, $method);
my @All_methods = qw(mean median std range min max);
my $matched;
foreach my $method_cur (@method){
	$matched = check_paras(\@All_methods, $method_cur);
	die "Unknown method: $method_cur!\n" unless ($matched == 1);
}

########################################################################################## 
# Output files
##########################################################################################
# Make output directary 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $basename =  basename($in_file, qw(.txt));
my $name = '';
if ($basename =~ /DicedNSIDE_(.*)/){
	$name = "dice_$1";
}elsif ($basename =~ /ConNSIDE_(.*)/){
	$name = "consensus_$1";
}elsif ($basename =~ /NSIDE_(.*)/){
	$name = $1;
}

my $output_annotation_file = '';
my $tmp_NSIDE_bed_file = '';
my $tmp_intersected_bed_file = '';

if ($name ne ''){
  $output_annotation_file = $out_folder . "/NSIDE_annotations_" . $name . ".txt";
  $tmp_NSIDE_bed_file = $out_folder . "/tmp_NSIDE_" . $name . ".bed";
  $tmp_intersected_bed_file = $out_folder . "/tmp_Intersected_" . $name . ".bed";
}else{
  $output_annotation_file = $out_folder . "/NSIDE_annotations.txt";
  $tmp_NSIDE_bed_file = $out_folder . "/tmp_NSIDE.bed";
  $tmp_intersected_bed_file = $out_folder . "/tmp_Intersected.bed";
}

########################################################################################## 
# Generate tmp bed file of NSIDE
##########################################################################################
my %NSIDE;

open ( NSIDE, "<", $in_file ) or die "$!";
open ( TMPNSIDEBED, ">", $tmp_NSIDE_bed_file) or die "$!";

while(<NSIDE>){

	next if (/Size/);

	chomp;
	my @F = split /\s+/;
	
	my @BeadIDs = split(/,/, $F[3]);
	my $NSIDEid = $F[0];

	foreach my $BeadID (@BeadIDs){
		my @site = split(/:|-/, $BeadID);
		my @chr = split(/_/, $site[0]);
		$site[0] = $chr[0];
		my $site = join "\t", @site;
		print TMPNSIDEBED "$site\t$NSIDEid\n";
	}

	$NSIDE{$NSIDEid}{Data} = join "\t", @F[1..2];
}

close NSIDE;
close TMPNSIDEBED;

########################################################################################## 
# Get column names for "direct" mode
##########################################################################################
my $names_of_each_column;

if ( $mode eq "direct"){

	my $header;

	if ($header_bed){
		open (BED, "<", $bed_file) or die "$!";
		$header = readline BED;	
		close BED;
	}

	$names_of_each_column = name_columns_by_methods(\@column_indexs, \@method, $header);
}

########################################################################################## 
# Intersect NSIDE bed file with input bed file
##########################################################################################
my $intersect_cmd;

if ($header_bed){
	$intersect_cmd = "sed '1d' $bed_file | intersectBed -a $tmp_NSIDE_bed_file -b - -f $fraction -wo > $tmp_intersected_bed_file";
}else {
	$intersect_cmd = "intersectBed -a $tmp_NSIDE_bed_file -b $bed_file -f $fraction -wo > $tmp_intersected_bed_file";
}

run_command($intersect_cmd); 

########################################################################################## 
# Annotate each NSIDE with elements ("element" mode) or values ("direct" mode)
##########################################################################################
my %ValueData;

open (INTERSECT, "<", $tmp_intersected_bed_file) or die "$!";
while(<INTERSECT>){
	chomp;
	my @F = split /\s+/;

	my $NSIDEid = $F[3];

	if ($mode eq "element"){
		my $E_name = $F[7];

		if ($NSIDE{$NSIDEid}{ElementsArray}){
			push @{$NSIDE{$NSIDEid}{ElementsArray}}, $E_name;
		}else{
			$NSIDE{$NSIDEid}{ElementsArray} = [($E_name)];
		}
	}elsif($mode eq "direct"){
		foreach my $column_index (@column_indexs){
			my $column_index_intersect = $column_index + 4;
			if ($ValueData{$NSIDEid}{$column_index}{ValueArray}){
				push @{$ValueData{$NSIDEid}{$column_index}{ValueArray}}, $F[$column_index_intersect];
			}else{
				$ValueData{$NSIDEid}{$column_index}{ValueArray} = [($F[$column_index_intersect])];
			}
		}
	}
}
close INTERSECT;

########################################################################################## 
# If input annotation file, record values of each element
##########################################################################################
my %ElementData;

if ($annotate_file){
	open (DATA, "<", $annotate_file) or die "$!";

	my $header;
	if ($header_annotate){
		$header = readline DATA;
	}
	$names_of_each_column = name_columns_by_methods(\@column_indexs, \@method, $header);

	while(<DATA>){
		chomp;
		next if /^\s*$/;
		my @F = split /\s+/;

		my $Element_id = $F[0];

		foreach my $column_index (@column_indexs){
			$ElementData{$Element_id}{$column_index} = $F[$column_index];
		}
	}
	close DATA;
}


########################################################################################## 
# Calculate values acoording to methods and output the results
##########################################################################################
open (OUT, ">", $output_annotation_file) or die "$!";

# Header of output file
if ($mode eq "element"){
	if ($detail){
		if ($names_of_each_column){
			print OUT "NSIDEid\tSize\tType\tElementNumber\t$names_of_each_column\tElementList\n";
		}else{
			print OUT "NSIDEid\tSize\tType\tElementNumber\tElementList\n";
		}
	}else{
		if ($names_of_each_column){
			print OUT "NSIDEid\tSize\tType\tElementNumber\t$names_of_each_column\n";
		}else{
			print OUT "NSIDEid\tSize\tType\tElementNumber\n";
		}	
	}
}elsif ($mode eq "direct"){
	print OUT "NSIDEid\tSize\tType\t$names_of_each_column\n";
}
	
# Calculate and output value of each NSIDE
foreach my $NSIDEid (sort keys %NSIDE){
	if ($mode eq "element"){
		my @Elements = @{$NSIDE{$NSIDEid}{ElementsArray}};
		@Elements = uniq(@Elements);
		my $Elements_num = @Elements;

		# Fill %ValueData
		foreach my $column_index (@column_indexs){
			my @values;
			foreach my $Element_id (@Elements){
				if ($ElementData{$Element_id}{$column_index}){
					push @values, $ElementData{$Element_id}{$column_index};
				}else{
					push @values, 0;
				}
			}
			$ValueData{$NSIDEid}{$column_index}{ValueArray} = [@values];
		}

		# Calculate
		my $values_of_each_column = calculate_columns_by_methods(\@column_indexs, \@method, \%ValueData, $NSIDEid);

		# Output
		if ($detail){
			my $ElementList = join ",", @Elements;

			if ($names_of_each_column){
				print OUT "$NSIDEid\t$NSIDE{$NSIDEid}{Data}\t$Elements_num\t$values_of_each_column\t$ElementList\n";
			}else{
				print OUT "$NSIDEid\t$NSIDE{$NSIDEid}{Data}\t$Elements_num\t$ElementList\n";
			}
		}else{
			if ($names_of_each_column){
				print OUT "$NSIDEid\t$NSIDE{$NSIDEid}{Data}\t$Elements_num\t$values_of_each_column\n";
			}else{
				print OUT "$NSIDEid\t$NSIDE{$NSIDEid}{Data}\t$Elements_num\n";
			}	
		}
	}elsif ($mode eq "direct"){
		#Calculate
		my $values_of_each_column = calculate_columns_by_methods(\@column_indexs, \@method, \%ValueData, $NSIDEid);

		# Output
		print OUT "$NSIDEid\t$NSIDE{$NSIDEid}{Data}\t$values_of_each_column\n";
	}		
}

##########################################################################################
# Clean up tmp file
##########################################################################################
unlink $tmp_NSIDE_bed_file;
unlink $tmp_intersected_bed_file;

##########################################################################################
# Subroutines
##########################################################################################
# Get 0-based columns from 1-based $column_str
# Parameters：$column_str
sub get_columns {
	my $column_str = shift; #1-based

	my @columns = ();   #0-based
	my %columns_tmp = ();

	my @column_strs = split /,/, $column_str;
	foreach my $column (@column_strs) {
		my ($start, $end ) = split /:/, $column ;
		$end = $start unless $end;
		my @tmp = ();
		my $revere_flag = 0;
		if ($start > $end) {
			($start, $end) = ($end, $start);
			$revere_flag = 1;
		}
		foreach my $col ($start .. $end) {
			if (exists $columns_tmp{$col}) {
				next;
			}else{
				$columns_tmp{$col} = 1;
				push @tmp, $col - 1;
			}
		}
		@tmp = reverse @tmp if $revere_flag;
		push @columns, @tmp;			
	}

	return @columns;
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

# Name columns by methods
# Parameters: \@column_indexes, \@methods, $header (optional)
# Notice: header is expected to be separated by tab
sub name_columns_by_methods {
	my @column_indexs =  @{$_[0]};
	my @methods = @{$_[1]};
	my $header = $_[2];

	my @names_of_each_column;

	if ($header){
		$header =~ s/\r?\n//;
		my @column_names = split(/\t/, $header); 

		foreach my $column_index (@column_indexs){
			foreach my $method (@methods){
				my $column_name = $column_names[$column_index];
				push @names_of_each_column, "$column_name\_$method";				
			}
		}
	}else{
		foreach my $column_index (@column_indexs){
			foreach my $method (@methods){
				my $column_name = $column_index + 1;
				push @names_of_each_column, "Col$column_name\_$method";
			}
		}	
	}
	my $names_of_each_column = join "\t", @names_of_each_column;
	return $names_of_each_column;
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


# For each NSIDE, calculate value in sepecified columns and methods
# Parameters: \@column_index, \@methods, \%ValueData, $NSIDEid
# Notice: %ValueData is expeted to be $ValueData{$NSIDEid}{$column_index}{ValueArray}
sub calculate_columns_by_methods {
	my @column_indexs =  @{$_[0]};
	my @methods = @{$_[1]};
	my %ValueData = %{$_[2]};
	my $NSIDEid = $_[3];

	my @results;
	foreach my $column_index (@column_indexs){
		my @values = @{$ValueData{$NSIDEid}{$column_index}{ValueArray}};

		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@values);

		foreach my $method (@methods){
			my $result_tmp;

			if ($method eq "mean"){
				$result_tmp = $stat->mean();
			}elsif ($method eq "median"){
				$result_tmp = $stat->median();
			}elsif ($method eq "std"){
				$result_tmp = $stat->standard_deviation();
			}elsif ($method eq "range"){
				$result_tmp = $stat->sample_range();
			}elsif ($method eq "min"){
				$result_tmp = $stat->min();
			}elsif ($method eq "max"){
				$result_tmp = $stat->max();
			}
			push @results, $result_tmp;
		}
	}

	my $values_of_each_column = join "\t", @results;
	return $values_of_each_column;
}
