#!/usr/bin/perl -w
use strict;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Statistics::R;
use Array::Utils;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_SemSim.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to quantitatively compute average
gene simularity (semantic simularity in directed
acyclic graph of gene ontology terms) inside single
NSIDE by GOSemSim and test whether the simularity is
higher than random in each slice of NSIDE_annotation 
file filt by the size of NSIDE.

Reference:
1) Yu et al., Methods in Molecular Biology (2020)
2) Yu et al., Bioinformatics (2010)

------- ------- ------- ------- ------- ------- -------

Iutput files (details in "Options"):

1. FiltInputSliceN generated by NSIDE_SemSim.pl,
   annotated NSIDE, including at least 5 columns:
   1)   NSIDEid;
   2)   Size of NSIDE: number of beads in the NSIDE;
   3)   Type of NSIDE: "Cis" means all beads in the
        NSIDE belong to same chromatin, otherwise
        "Trans";
   4)   List of BeadID of all beads in the NSIDE.
   5)   Number of elements;
   6-n) If input annotation file: specified
        statistics of values (by --method) in each
        specfied column (by -c).
   n+1) List of elements.
2. GenesAvailable_ont_drop.ids: id of all annotated
   genes in the selected aspect of GO after dropping
   the selected evidence code(s). 

------- ------- ------- ------- ------- ------- -------

Output file:

1. OutSemN: annotated NSIDE containing semantice
   simularity results:
   1-n) Same as FiltInputSliceN;
   n+1) SemSim: average semantic simularity of genes;
   n+2) zscore_SemSim: standardized z of semantic
        simularity; 
   n+3) p_SemSim: pvalue of semantic simularity;
   > If "--readable":
   n+4) Official gene symbols of genes;
   n+5) List of elements.
   > Else:
   n+4) List of elements.

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_SemSim_sub.pl -i FiltInputSliceN
                    --gene GenesAvailable_ont_drop.ids
                    [Options]

NSIDE_SemSim_sub.pl -i FiltInputSliceN
                    --gene GenesAvailable_ont_drop.ids
                    --mean 0.3913242 --sd 0.5449009
                    --species Mm --keytype REFSEQ
                    --ont BP --drop NULL 
                    --measure Wang --combine max
                    --readable
                    -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:
 -i|--input      > FiltInputSliceN file.
                 > !Notice:
                 > Must contain a list of genes in the
                   last column!
                 > The keytype of gene names should NOT
                   be "SYMBOL", "ALIAS", "COMMON" or
                   "GENENAME" because gene symbols
                   cannot be relied upon to be uniquely
                   mapped onto a single gene.
 --gene          > GeneAvailable.ids file.
 --mean          > Mean simularity of all gene pairs,
                   estimated by NSIDE_SemSim.pl.
 --sd            > Standard deviation of all gene
                   pairs, estimated by NSIDE_SemSim.pl.
 --species       > Species, "Mm" for house mouse
                   (default) or "Hs" for human.
 --keytype       > Keytype of gene names. "REFSEQ"
                   (default), "ACCNUM", "ENSEMBL",
                   "ENSEMBLPROT", "ENSEMBLTRANS",
                   "ENTREZID", "ENZYME", "INTERPRO", 
                   "IPI", "OMIM", "PROBEID", "UNIGENE",
                   or "UNIPROT", please refer to the
                   mananual of R package
                   "AnnotationDbi".
                 > !Notice:
                 > The keytype of gene names should NOT
                   be "SYMBOL", "ALIAS", "COMMON" or
                   "GENENAME".
 --ont           > Aspect of ontology, "BP" for
                   biological process (default), "CC"
                   for cellular component or "MF" for
                   molecular function.
 --drop          > One or a set of evidence codes
                   separated by "," to drop from
                   calculating semantic simularity.
                   Default is NULL, please refer to
                   http://geneontology.org/docs/
                   guide-go-evidence-codes/
 --measure       > One of "Resnik", "Lin", "Rel",
                   "Jiang" and "Wang" (default).
                 > Notice!:
                 > If any evidence codes are dropped,
                   use IC-based methods (including
                   "Resnik", "Lin" and "Rel") with
                   caution because the IC value may
                   not be accurate!
                 > Please refer to https://yulab-smu.
                   top/biomedical-knowledge-mining-
                   book/semantic-similarity-overview.
                   html
 --combine       > One of "max" (default), "avg",
                   "rcmax", "BMA" for combining
                   simularity scores of multiple GO
                   tems associated with one gene.
                 > Notice!:
                 > "BMA" may lead to error in some
                   situations, if so, please try
                   "max", "avg" or "rcmax".
                 > If use "--ont CC", "avg" is better
                   than "max" in estimating mean and
                   sd  of simularity.
 --readable      > "TRUE" or "FALSE" (default), conver
                   the gene names in the last column
                   into official gene symbols and add a
                   column in the output file.
                 > !Notice:
                 > The number of official gene symbols
                   may not be equal to the
                   "ElementNumber" because multiple
                   genes may have same official gene
                   symbol.
 -o|--output     > Output folder.
 -h|--help       > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my $in_file    = '';
my $gene_file  = '';
my $mean       = 0.3913242;
my $sd         = 0.5449009;
my $species    = 'Mm';
my $keytype    = 'REFSEQ';
my $ont        = 'BP';
my $drop       = 'NULL';
my $measure    = 'Wang';
my $combine    = 'max';
my $random     = 100;
my $force      = 0;
my $readable   = "FALSE";
my $out_folder = dirname './';
my $help       = 0;

die $usage
  unless GetOptions(
  	"i|input=s"       => \$in_file,
  	"gene=s"          => \$gene_file,
  	"mean:f"          => \$mean,
  	"sd:f"            => \$sd,
    "species:s"       => \$species,
  	"keytype:s"       => \$keytype,
  	"ont:s"           => \$ont,
  	"drop:s"          => \$drop,
  	"measure:s"       => \$measure,
  	"combine:s"       => \$combine,
  	"random:i"        => \$random,
   	"force"           => \$force,
  	"readable:s"      => \$readable,
    "o|output:s"      => \$out_folder,
    "h|help"          => \$help,
  );

##########################################################################################
# Check the parameter infomation 
##########################################################################################  
die $usage if $help;

# Check the input file
die "No FiltInputSliceN file!\n" unless -s $in_file;
die "No GenesAvailable.ids file!\n" unless -s $gene_file;

# Check the parameters
# --mean and --sd
die "No mean!\n" unless $mean;
die "No sd!\n" unless $sd;

# --species
my $library;
my $bimap;
my $srcSpecies;
if ($species eq "Hs"){
	$library = 'org.Hs.eg.db';
	$bimap = 'org.Hs.egGO';
	$srcSpecies = 'HOMSA';
}elsif ($species eq "Mm"){
	$library = 'org.Mm.eg.db';
	$bimap = 'org.Mm.egGO';
	$srcSpecies = 'MUSMU';
}

# --drop
my @drop = split(/,/, $drop);
my $drop_R = join(",", map { qq/"$_"/ } @drop);

# --measure
my $computeIC;
if ($measure eq "Wang"){
	$computeIC = "FALSE";
}elsif ("Resnik|Lin|Rel|Jiang" =~ $measure){
	$computeIC = "TRUE";
}elsif ($measure eq "TCSS"){
	$computeIC = "TRUE";
}

# --readable
die "--readable is expected to be \"TRUE\" or \"FALSE\"!\n" unless ($readable eq "TRUE" | $readable eq "FALSE");

########################################################################################## 
# Output files
##########################################################################################
# Make output directary
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $basename =  basename($in_file);
my $name = '';
if ($basename =~ /FiltInputSlice(.*)/){$name = $1;}

my $out_semsim_file = $out_folder . '/OutSem' . $name;

##########################################################################################
# Load the DAG of GO and get available gene ids
##########################################################################################
my $datestring = localtime(); 
warn("$datestring: For $in_file, loading Map...\n");

my $R = Statistics::R->new();

my $R_cmds = <<EOF;
	suppressPackageStartupMessages(library(GOSemSim))
	suppressPackageStartupMessages(library($library))
	GOMap <- godata("$library", keytype = "$keytype", ont="$ont", computeIC = "$computeIC", processTCSS = "FALSE")
	mean_random <- $mean
	sd_random <- $sd
	mean_random
	sd_random
EOF

$R -> run($R_cmds);

##########################################################################################
# Build vocab of available genes
##########################################################################################
$datestring = localtime(); 
warn("$datestring: For $in_file, building vocab of available genes...\n");

my %Available;
open( GA, "<", $gene_file ) or die "$!";
while(<GA>){
	chomp;
	$Available{$_} = 1;
}
close GA;

##########################################################################################
# Calculate semantic simularity of each NSIDE
##########################################################################################
$datestring = localtime(); 
warn("$datestring: For $in_file, calculating semantic simularity of each NSIDE...\n");

open( OUTPUT, ">", $out_semsim_file ) or die "$!";
open( INPUT, "<", $in_file ) or die "$!";
while(<INPUT>){
	chomp;
	my @F = split /\t/;

	my $NSIDEid = $F[0];
	my $Size = $F[1];

	# $datestring = localtime(); 
	# warn("$datestring: For $in_file, processing $NSIDEid\...\n");

	my @Genes = split(/,/, $F[-1]);
	next if (@Genes <= 1);
	my $count_avai = 0;
	foreach my $gene (@Genes){
		if ($Available{$gene}){
			$count_avai += 1;
			if ($count_avai == 2){
				last;
			}
		}
	}

	if ($count_avai == 2){
		my $content = join "\t", @F[0..-2];
		my $Genes = join(",", map { qq/"$_"/ } @Genes);
		my $GeneNum = @Genes;
		my $n = ($GeneNum * ($GeneNum - 1)) * 0.5;

		my $R_cmds_cur = <<EOF;
		Genes = c($Genes)
		Mtx <-mgeneSim(Genes, semData=GOMap, drop = c($drop_R), measure="$measure", combine="$combine", verbose=FALSE)
		Result <- mean(Mtx[upper.tri(Mtx)])
		Zscore = (Result - mean_random)/(sd_random / sqrt($n)) 
		pvalue = pnorm(q=Zscore, lower.tail=FALSE)
		Result
		Zscore
		pvalue
		if ($readable){
			Genes_Readable <- intraIDMapper(Genes, species="$srcSpecies", srcIDType="$keytype", destIDType="SYMBOL")
			Genes_Readable <- as.vector(unlist(Genes_Readable))
			Genes_Readable
		}
EOF
		$R -> run($R_cmds_cur);
		my $Result = $R -> get("Result");
		my $Zscore = $R -> get("Zscore");
		my $pvalue = $R -> get("pvalue");
		my @content = ($Result, $Zscore, $pvalue);

		if ($readable eq "TRUE"){
			my $Genes_Readable_ref = $R -> get("Genes_Readable");
			if (ref $Genes_Readable_ref){
				my @Genes_Readable = @{$Genes_Readable_ref};
				@Genes_Readable = Array::Utils::unique(@Genes_Readable);
				my $GeneNum_Readable = @Genes_Readable;
				my $Genes_Readable = join ",", @Genes_Readable;
				push @content, $Genes_Readable;
			}else{
				push @content, $Genes_Readable_ref;
			}
		}

		splice(@F, -1, 0, @content);
		print OUTPUT join "\t", @F;
		print OUTPUT "\n";
	}
}
$R -> stop();
close INPUT;
close OUTPUT;

##########################################################################################
# Clean up tmp file
##########################################################################################
unlink $in_file;