#!/usr/bin/perl -w
use strict;
use 5.010;
use File::Basename;
use IPC::Cmd qw[can_run run];
use Array::Utils qw[array_diff];
use List::Util qw[sum];
use Statistics::Descriptive; 
use Statistics::R;
use Getopt::Long;

my $usage = <<USAGE;

SYSNOPSIS
------- ------- ------- ------- ------- ------- -------
NSIDE_LISA.pl V5.17, written by Xinyi Liu       
------- ------- ------- ------- ------- ------- -------

This program is used to calculate the local Moran’s I 
proposed by Anselin of input beads.

Reference:
Anselin Luc. Geographical Analysis (1995)

------- ------- ------- ------- ------- ------- -------

Input files:

1. RawData.txt: from NSIDE_Pre.pl, refomated and 
   annotated data of each bead (genomic bin), including 
   5 columns: BeadID, x, y, z, Value.
2. EdgeBeads.txt: from NSIDE_Edge.pl, edged beads of 
   each beads, including 2 columns: BeadID, list of 
   BeadID of all edged beads.

------- ------- ------- ------- ------- ------- -------

Output files:

1. LISAresults.txt: details of local Moran's I of each
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
2. SigBeads.txt: significant beads (FDR < cutoff, Type 
   = 1, columns are the same as LISAresults.txt).

------- ------- ------- ------- ------- ------- -------

Usage:

NSIDE_LISA.pl -i RawData_xx.txt -e EdgedBeads_xx.txt 
              [Options]

NSIDE_LISA.pl -i RawData_xx.txt -e EdgedBeads_xx.txt 
              -cf 0.05 -o output_folder_name

------- ------- ------- ------- ------- ------- -------

Options:

 -i|--inputfile  > RawData file from NSIDE_Pre.pl.
 -e|--edgefile   > EdgedBeads file from NSIDE_Edge.pl.
 -cf|--cutofffdr > Cutoff of FDR, default 0.05.
 -o|--output     > Output folder.
 -h|--help       > Print this usage.

------- ------- ------- ------- ------- ------- -------

USAGE
my $in_file         = '';
my $edge_file       = '';
my $cutofffdr       = 0.05;
my $out_folder      = dirname './';
my $help            = 0;

die $usage
  unless GetOptions(
    "i|inputfile=s"   => \$in_file,
    "e|edgefile=s"    => \$edge_file,
    "cf|cutofffdr:f"  => \$cutofffdr,
    "o|output:s"      => \$out_folder,
    "h|help"          => \$help,
  );

##########################################################################################
# Check the parameter infomation 
##########################################################################################  
die $usage if $help;

# Check the running environment
can_run('R') or die 'R is not executable!';

# Check the input file
die "No RawData file!\n" unless -s $in_file;
die "No EdgeBeads file!\n" unless -s $edge_file;

#check the cutoff value
die "Cutoff of FDR is expected to (0,1]!" unless ($cutofffdr > 0 && $cutofffdr <=1);

########################################################################################## 
# Output files
##########################################################################################
# Make output directary 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;

# File handles
my $basename =  basename($in_file, qw(.txt));
my $name = '';
if ($basename =~ /RawData_(.*)/){$name = $1;}

my $output_lisa = '';
my $output_sig_beads = '';

if ($name ne ''){
  $output_lisa = $out_folder . "/LISAresults_" . $name . ".txt";
  $output_sig_beads = $out_folder . "/SigBeads_" . $name . ".txt"; 
}else{
  $output_lisa = $out_folder . "/LISAresults.txt";
  $output_sig_beads = $out_folder . "/SigBeads.txt"; 
}

##########################################################################################
# Read RawData and count N (total beads number) and n (empty beads number)
# Get Xmean and Xsd
########################################################################################## 
my %BeadData = ();
my $Xmean;
my $Xsd;
my $N = 0;
my $n = 0;
my @Values = ();

open( INPUT, "<", $in_file ) or die "$!";
while(<INPUT>){
  next if (/BeadID/);

  if (/\A###/){
    if ($_ =~ /Xmean = (.*), Xsd = (.*)\Z/){ $Xmean = $1, $Xsd = $2};
  }else{
    chomp;

    my @F = split /\s+/, $_;

    my $BeadID = $F[0];
    my $Value = $F[4];

    if ($Value != 0){
      $BeadData{$BeadID}{Value} = $Value;
      push @Values, $Value;
    }else{
      $n += 1;
    }

    $N += 1;
  }
}
close INPUT;

##########################################################################################
# Calculate Xi (Normalized Value), Si2, b2i
##########################################################################################
my %ValueData;

foreach my $BeadID (keys %BeadData){

  my $Value = $BeadData{$BeadID}{Value};

  # Calculate Xi (Normalized Value)
  $BeadData{$BeadID}{Xi} = ($Value - $Xmean)/$Xsd;

  # Calculate Si2 and b2i  
  if ($ValueData{$Value}){
    $BeadData{$BeadID}{Si2} = $ValueData{$Value}{Si2};
    $BeadData{$BeadID}{b2i} = $ValueData{$Value}{b2i};  
  }else{
    my @value_tmp = ();
    push @value_tmp, $Value;
    my @allothervalue = array_diff(@Values, @value_tmp);

    my @numerator_forSi2 = ();
    my @numerator_forb2i = ();
    my @denominator_forb2i = ();

    foreach (@allothervalue){
        push @numerator_forSi2, ($_ - $Xmean)**2;
        push @numerator_forb2i, ($_ - $Xmean)**4;
        push @denominator_forb2i, ($_ - $Xmean)**2;
    }

    my $Si2 = ((sum @numerator_forSi2) + ($n* ($Xmean)**2))/($N - 1);
    my $b2i = ((sum @numerator_forb2i) + ($n* ($Xmean)**4))/(((sum @denominator_forb2i) + ($n* ($Xmean)**2))**2);

    $BeadData{$BeadID}{Si2} = $Si2;
    $BeadData{$BeadID}{b2i} = $b2i;

    $ValueData{$Value}{Si2} = $Si2;
    $ValueData{$Value}{b2i} = $b2i; 
  }
}

##########################################################################################
# Read edge infomation and calculate Ii, wXj A, B, EI2, VI, zI, LISA, Cluster
########################################################################################## 
open( EDGE, "<", $edge_file) or die "$!";
open( LISA, ">", $output_lisa ) or die "$!";
print LISA "BeadID\tValue\tIi\tzI\tXi\twXj\tCluster\n";

while(<EDGE>){
  next if (/BeadID/);
  chomp;
  
  my @F = split /\s+/;
  my $BeadID = $F[0];

  next unless ($BeadData{$BeadID}{Value});

  if ($F[1]){
    my @EdgeBeadIDs = split(/,/, $F[1]);

    # Calculating Ii and wXj (sXj (Spatial Lag))
    my @sthneedsum_forIi = ();
    my @EdgeCubewXj = ();

    foreach my $EdgeID (@EdgeBeadIDs){ 
      if (defined $BeadData{$EdgeID}{Value}){
          push @sthneedsum_forIi, ($BeadData{$EdgeID}{Value} - $Xmean);
          push @EdgeCubewXj, $BeadData{$EdgeID}{Xi};
      }else{
          push @sthneedsum_forIi, (0 - $Xmean);
          push @EdgeCubewXj, 0;
      }
    }
    
    # Ii  
    my $Ii = ((sum @sthneedsum_forIi)*($BeadData{$BeadID}{Value} - $Xmean))/$BeadData{$BeadID}{Si2};

    # wXj
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@EdgeCubewXj);
    my $wXj = $stat->mean();

    # Calculating A, B and E[I**2]（EI2）
    my $a = @EdgeBeadIDs; #if edged, wij=1, sum(wij^2)
    my $b = (@EdgeBeadIDs)**2; #sum(wik x wih (k!=i, h!=i))
    my $A = (($N - $BeadData{$BeadID}{b2i})*$a)/($N - 1);
    my $B = (((2*$BeadData{$BeadID}{b2i})-$N)*$b)/(($N- 1)*($N - 2));  
    my $EI2= $A + $B;

    # Calculating VI and zI
    my $EI = 0 - (@EdgeBeadIDs / ($N - 1));
    my $VI = $EI2 - ($EI**2);
    my $zI = ($Ii - $EI)/sqrt($VI);

    # Type according to Xi and wXj (sXj)
    # HHHH（++）==1 HLLH（-+）==2
    # LLLL（--）==3 LHHL（+-）==4
    my $Type;
    if ($BeadData{$BeadID}{Xi} > 0 && $wXj > 0){
        $Type = 1;
    }elsif($BeadData{$BeadID}{Xi} < 0 && $wXj < 0){
        $Type = 3;
    }elsif($BeadData{$BeadID}{Xi} < 0 && $wXj > 0){
        $Type = 2;
    }elsif($BeadData{$BeadID}{Xi} > 0 && $wXj < 0){
        $Type = 4;   
    }else{
        $Type = 0;
    }
    
    print LISA "$BeadID\t$BeadData{$BeadID}{Value}\t$Ii\t$zI\t$BeadData{$BeadID}{Xi}\t$wXj\t$Type\n";
    
  }else{
    # For beads without edged beads:
    # Ii = 0 (na), zI = 0 (na), wXj = 0 (na), Type == 5 
    print LISA "$BeadID\t$BeadData{$BeadID}{Value}\t0\t0\t$BeadData{$BeadID}{Xi}\t0\t5\n";    
  }
}
close EDGE;

########################################################################################## 
# Calculate p value and FDR
##########################################################################################
my $R = Statistics::R->new();
my $cmds_R = <<EOF;
l.moran <- data.frame(read.table(file = "$output_lisa", header = TRUE))
l.moran[,8] <- 2*pnorm(-abs(l.moran[,4]))
l.moran[,9] <- p.adjust(l.moran[,8], method="bonferroni")
l.moran[,10] <- p.adjust(l.moran[,8], method="fdr")
l.moran
names(l.moran) <- c("BeadID", "Value", "Ii", "zI", "Xi", "wXj", "Cluster", "pvalue", "p_bonferroni", "FDR")
write.table(l.moran, file = "$output_lisa" , quote = FALSE, append = FALSE, sep = "\t", na = "NA", dec = ".", row.name = FALSE, col.names = TRUE, fileEncoding = "")
SigBeads <- subset(l.moran, FDR < 0.05 & Cluster == 1)
write.table(SigBeads, file = "$output_sig_beads" , quote = FALSE, append = FALSE, sep = "\t", na = "NA", dec = ".", row.name = FALSE, col.names = TRUE, fileEncoding = "")
EOF

$R->run($cmds_R);
