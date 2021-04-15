#!/usr/bin/env perl

use strict;
#use warnings;
use Getopt::Std;
use List::Util qw(shuffle);

my $PRGM_NAME = 'quant-noise.pl';
my $VERSION_DATE = '10/18/17';
my $SIMBASE = 'sim_noise';
my $CHANGEFRACT = 0.3;


# replace the second SNP in the pedigree file with a set of noise SNPs
# create a trait file with the correlated phenotype and then the noise phenotypes

sub usage{
   print "\n\t$PRGM_NAME:  $VERSION_DATE\n\n";
   print "\tusage:\t\t$PRGM_NAME -f <ped file> -m <MAF> -s <num SNPs> -p <num Phenotypes> -r <random seed>\n\n";
   print "\texample:\t$PRGM_NAME -f input.ped -m 0.1 -s 3 -p 9 -r 1799\n\n";

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
   my $pedfile = $opts{f} || die "Requires pedigree filename (-f)\n\n";
   my $randseed = $opts{r} || 777;

   print "num Phenotypes=$nPhen\n";
   print "num SNPs=$nSnps\n";
   print "maf=$maf\n";
   print "pedfile=$pedfile\n";
   print "randseed=$randseed\n";

   return ($pedfile, $maf, $nSnps, $nPhen, $randseed);
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


# add noise SNPs with MAF specified
sub create_newped{
  my $inds = shift;
  my $maf = shift;
  my $nSNPs = shift;
  my $nPhenos=shift;

#     add the Phneotype noise created by changing a fraction (0.3) of the original genotypes
#     multiply the num of inds by the fraction and then choose that number to mutate
#     change to one of other genotypes based on odds 
    my @noisePhenos;
    for(my $ph=0; $ph<$nPhenos; $ph++){
      push(@{$noisePhenos[$ph]},create_pheno_noise2($inds,$maf));
    }

  my $newped = "$SIMBASE.ped";
  open(OUT, ">$newped") or die "$newped:  $!\n\n";
  print OUT "#IID Pheno SNP1_1 SNP1_2";
  for(my $i=2; $i<=$nSNPs+$nPhenos+1; $i++){
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

#    for(my $p=0; $p<@noisePhenos; $p++){
    foreach my $phenoarray(@noisePhenos){
       print OUT " $$phenoarray[$i][0] $$phenoarray[$i][1]";      
    }

    print OUT "\n";
  }

  close(OUT); 
}

# create phenotype noise by changing a fraction of the original phenotypes
sub create_pheno_noise1{
  my $inds = shift;
  my $maf = shift;
  my $nChange = int($CHANGEFRACT*scalar(@$inds));

  my @genos; 
  $genos[1][0]=2;
  $genos[1][1]=2;
  $genos[2][0]=1;
  $genos[2][1]=1;
  $genos[3][0]=1;
  $genos[3][1]=2;
  $genos[4][0]=2;
  $genos[4][1]=2;
  $genos[5][0]=1;
  $genos[5][0]=1;
  
  my @addchance;
  $addchance[2] = $maf*(1-$maf)*2 / ($maf*(1-$maf)*2 + $maf*$maf);
  $addchance[3] = $maf*$maf / ($maf*$maf + ((1-$maf)*(1-$maf)));
  $addchance[4] = (1-$maf)*(1-$maf) / ($maf*(1-$maf)*2 + ((1-$maf)*(1-$maf)));

  my @indgenos;
  my @indexes;
  for(my $i=0; $i<@$inds; $i++){
    $indgenos[$i] = [$$inds[$i][2], $$inds[$i][3]];
    push(@indexes,$i);
  }
  
  my @shuffled = shuffle(@indexes);
  my %selected;
  for(my $i=0; $i<$nChange; $i++){
    $selected{$shuffled[$i]}=1;
  } 

  # when constructing the array for change only change those found in the selected hash
  my @newgenos;
  for(my $i=0; $i<@indgenos; $i++){
    if($selected{$i}){
      my $geno;
      $geno=$indgenos[$i][0]+$indgenos[$i][1];
      if(rand() < $addchance[$geno]){
        push(@newgenos, $genos[$geno+1]);
      }
      else{
        push(@newgenos, $genos[$geno-1]);
      }
    }
    else{
      push(@newgenos, $indgenos[$i]);
    }
  }
  return @newgenos;  
}


# alternate method for changing genotypes portion (0.3)
# swap values of 0.3 samples -- MAF is unchanged in this case
sub create_pheno_noise2{
  my $inds = shift;
  my $maf = shift;
  my $nChange = int($CHANGEFRACT*scalar(@$inds));

  my @indexes;
  my @indgenos;
  for(my $i=0; $i<@$inds; $i++){
    $indgenos[$i] = [$$inds[$i][2], $$inds[$i][3]];
    push(@indexes,$i);
  }
  my @shuffled = shuffle(@indexes);
  # swap genotypes with next individual in shuffled array
  # last is swapped for first
  my $firstgenos=$indgenos[$shuffled[0]];
  my $idx;
  for($idx=0; $idx<$nChange-1; $idx++){
    $indgenos[$shuffled[$idx]]=$indgenos[$shuffled[$idx+1]]; 
#print "switch $shuffled[$idx] with $shuffled[$idx+1]\n";
  }
  $indgenos[$shuffled[$idx]]=$firstgenos;
#print "switch $shuffled[$idx] with $shuffled[0]\n";
#print "-------\n";
  return @indgenos;
}

sub create_map{
  my $nSNPs=shift;
  my $nPhenos=shift;
  my $newmap="$SIMBASE.map";
  open(OUT, ">$newmap") or die "$newmap:  $!\n\n";
  print OUT "1 rsCorr 1 1 2\n";
  my $i;
  for($i=1; $i<=$nSNPs; $i++){
    print OUT "1 rsNoise$i ", $i+1, " 1 2\n";
  }
  for(my $j=1; $j<=$nPhenos; $j++){
    my $num=$i+$j;    
    print OUT "1 rsPhenoNoise$num ",$i+$j, " 1 2\n";
  }
  close(OUT);
}


my ($pedfile, $maf, $nSNPS, $nPhenos, $randseed) = check_args();
srand($randseed);
my $inds=parse_pedfile($pedfile);
create_newped($inds, $maf, $nSNPS,$nPhenos);
create_map($nSNPS, $nPhenos);

