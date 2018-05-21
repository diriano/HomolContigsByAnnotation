#!/usr/bin/perl

use strict;
use warnings;
use Bio::Range;

###################################
## This script will read LASTZ files
##  produced with the option
##  --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,nmatch,nmismatch,ngap,continuity,identity,coverage
###################################
die if !@ARGV;

my $minHSPcov=1;
my $minHSPident=50;
my @ranges=();
my %rangeSeen;

foreach my $file(@ARGV){
 open INFILE, $file;
 my $largestContig='';
 my $smallestContig='';
 while(<INFILE>){ 
  chomp;
  next if /^#/;
  my ($score,$name1,$strand1,$size1,$zstart1,$end1,$name2,$strand2,$size2,$zstart2,$end2,$nmatch,$nmismatch,$ngap,$continuity,$conPct,$identity,$idPct,$coverage,$covPct)=split(/\t/);
  $largestContig=$name1;
  $smallestContig=$name2;
  if ($size2 > $size1){ die "The second contigs should be the smaller one.\n";}
  $idPct=~s/%//;
  $covPct=~s/%//;
  if($covPct >= $minHSPcov && $idPct >= $minHSPident){
#   print $_."\n";
   push @ranges, Bio::Range->new(-start=>$zstart2, -end=>$end2, -strand=>$strand2);
  }
 }
 my $nonOverlapingRanges=reduceRanges(\@ranges,1);
 foreach my $r(@$nonOverlapingRanges){
  print "$smallestContig ".$r->toString()."\n";
  ###TODO: compute aligned fraction, and take desicion whether or not to send to haplotig pool
 }
}

sub reduceRanges{
 my $ran=shift;
 my @ran=@$ran;
 my $overlapCount=shift;
 my @ranNew=();
 my $overlapsPresent=0;
# print "XXX ".scalar(@ran)." $overlapCount\n\n";
 if($overlapCount>0){
  for (my $i=0;$i<@ran;$i++){
#   print "i=$i\n\n";
   my $overlapCountNew=0;
   next if $rangeSeen{$ran[$i]->start}{$ran[$i]->end};
   for (my $j=$i+1;$j<@ran;$j++){
    next if $rangeSeen{$ran[$j]->start}{$ran[$j]->end};
    if($ran[$i]->overlaps($ran[$j])){
     $overlapCountNew++;
     $overlapsPresent++;
     my ($start, $stop, $strand) = $ran[$i]->union($ran[$j]);
     push @ranNew, Bio::Range->new(-start=>$start, -end=>$stop, -strand=>$strand);
# print "YYYXXX $overlapCount $overlapCountNew\n\n";
     $rangeSeen{$ran[$j]->start}{$ran[$j]->end}=1;
#     print $ran[$i]->start ."\t".$ran[$i]->end."\n";
#     print "\t".$ran[$j]->start ."\t".$ran[$j]->end."\n";
#     print "\tNEWRANGE= $start, $stop, $strand\n";
    }
   }
   if($overlapCountNew==0){
#    print "DIEGO\n";
    push @ranNew, Bio::Range->new(-start=>$ran[$i]->start, -end=>$ran[$i]->end, -strand=>$ran[$i]->strand);
   }
  }
#  print "RANGES IN NEW= ".scalar(@ranNew)."\n";
  reduceRanges(\@ranNew,$overlapsPresent);
 }
 else{
#  print "FINAL RANGES IN NEW= ".scalar(@ran)."\n";
#  foreach my $r(@ran){
#   print "PP ".$r->toString()."\n";
#  }
#  print "FINISHED\n\n";
  return \@ran
 }
}
