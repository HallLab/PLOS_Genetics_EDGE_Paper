#!/usr/bin/env perl

use strict;
#use warnings;
use Getopt::Std;
use List::Util qw(shuffle);

my $PRGM_NAME = 'casecon-noise.pl';
my $VERSION_DATE = '10/11/17';
my $SIMBASE = 'sim_noise';


# replace the second SNP in the pedigree file with a set of noise SNPs
# create a trait file with the correlated phenotype and then the noise phenotypes

sub usage{
   print "\n\t$PRGM_NAME:  $VERSION_DATE\n\n";
   print "\tusage:\t\t$PRGM_NAME -f <ped file> -m <MAF> -s <num SNPs> -c <Cases> -p <num Phenotypes> -r <random seed>\n\n";
   print "\texample:\t$PRGM_NAME -f input.ped -m 0.1 -s 3 -c 100 -p 9 -r 1799\n\n";

}

sub check_args{
   if(@ARGV < 4){
      usage();
      exit;
   }

   my %opts;
   getopts('f:m:s:c:p:r:',\%opts);

   my $nPhen = $opts{p} || 9;
   my $nSnps = $opts{s} || 3;
   my $maf = $opts{m} || die "Requires maf value (-m)\n\n";
   my $nCases = $opts{c} || die "Requires number cases (-c)\n\n";
   my $pedfile = $opts{f} || die "Requires pedigree filename (-f)\n\n";
   my $randseed = $opts{r} || 777;

   print "num Phenotypes=$nPhen\n";
   print "num SNPs=$nSnps\n";
   print "maf=$maf\n";
   print "num Cases=$nCases\n";
   print "pedfile=$pedfile\n";
   print "randseed=$randseed\n";

   return ($pedfile, $maf, $nSnps, $nCases, $nPhen, $randseed);
}


sub parse_pedfile{
  my $filename = shift;
  open(IN, $filename) or die "$filename:  $!\n\n";

  my @inds;
  my $headerline = <IN>;
  while(<IN>){
    chomp;
    my @info = (split" ")[0..3];
    push(@inds, \@info); 
  }
  close(IN);
  print "num of inds = ", scalar(@inds);
  print "\n";
  return \@inds;
}

sub create_traitfile{
  my $inds = shift;
  my $nCases = shift;
  my $nPhenos = shift;

  my @corr;
  for(my $i=0; $i<@$inds; $i++){
    push(@corr, $$inds[$i][1]);
  }

  my @phenos;
  push(@phenos, \@corr);

  for(my $i=1; $i<=$nPhenos; $i++){
    my @noise = shuffle(@corr); 
    push(@phenos, \@noise);
  }

  open(OUT, ">${SIMBASE}_traits.txt") or die "${SIMBASE}_traits.txt:  $!\n\n";
  # print header
  print OUT "IID CORR";
  for(my $i=1; $i<=$nPhenos; $i++){
    print OUT " NOISE$i";
  }
  print OUT "\n";

  my $nInds = scalar(@$inds);
  for(my $i=0; $i<$nInds; $i++){
    print OUT $i+1," $phenos[0][$i]";
    for(my $j=1; $j<=$nPhenos; $j++){
      print OUT " $phenos[$j][$i]";
    }
    print OUT "\n";
  }
  close(OUT);
}

# add noise SNPs with MAF specified
sub create_newped{
  my $inds = shift;
  my $maf = shift;
  my $nSNPs = shift;

  my $newped = "$SIMBASE.ped";
  open(OUT, ">$newped") or die "$newped:  $!\n\n";
  print OUT "#IID Pheno SNP1_1 SNP1_2";
  for(my $i=2; $i<=$nSNPs+1; $i++){
    print OUT " SNP${i}_1 SNP${i}_2";
  }
  print OUT "\n";
  for(my $i=0; $i<@$inds; $i++){
    printf OUT "%d $$inds[$i]->[1] $$inds[$i]->[2] $$inds[$i]->[3]", $$inds[$i]->[0];
    for(my $j=0; $j<$nSNPs; $j++){
      if(rand() < $maf){
        if(rand() < $maf){
          print OUT " 2 2";
        }
        else{
          print OUT " 1 2";
        }
      } 
      else{
        if(rand() < $maf){
          print OUT " 1 2";
        } else{
          print OUT " 1 1";
        }
      }
    }
    print OUT "\n";
  }

  close(OUT); 
}

sub create_map{
  my $nSNPs=shift;
  my $newmap="$SIMBASE.map";
  open(OUT, ">$newmap") or die "$newmap:  $!\n\n";
  print OUT "1 rsCorr 1 1 2\n";
  for(my $i=1; $i<=$nSNPs; $i++){
    print OUT "1 rsNoise$i ", $i+1, " 1 2\n";
  }
  close(OUT);
}


my ($pedfile, $maf, $nSNPS, $nCases, $nPhenos, $randseed) = check_args();
srand($randseed);
my $inds=parse_pedfile($pedfile);
create_traitfile($inds, $nCases, $nPhenos);
create_newped($inds, $maf, $nSNPS);
create_map($nSNPS);

