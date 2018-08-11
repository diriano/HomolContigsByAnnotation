#!/usr/bin/perl

use strict;
use warnings;
use Bio::Range;

###################################
## This script will read LASTZ files
##  produced with the option
##  --format=general:score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,nmatch,nmismatch,ngap,continuity,identity,coverage
##
## This script will filter the LASTZ alignments, first keepoing only High Scoring Pairs (HSP) with at least ($minHSPcov coverage or $minHSPAlgnBases aligned bases)  and at least $minHSPident identity
##  Then, it will reduce the range to merge all overlaping alignemnts and produce a final set of non-overlaping alingments, in order to be considered an
##  alternative hoplotype the final aligned length should be at least $minTotalFracAlgn % of the total length of the shortest contig
###################################
die if !@ARGV;

my $minHSPcov=3;
my $minHSPident=90;
my $minHSPAlgnBases=250;
my $minTotalFracAlgn=80;
my %contigSet1;
my %contigSet2;
my $fileCandidateContigsForRemoval="candidateContigs4Removal.txt";
my $fileContigs2HaplotigPool="Contigs2HaplotigPool.txt";

open CC4R, ">$fileCandidateContigsForRemoval";
open C2H, ">$fileContigs2HaplotigPool";


foreach my $file(@ARGV){
 open INFILE, $file;
 my $largestContig='';
 my $smallestContig='';
 my @ranges=();
 my %rangeSeen;
 my %lengths;

 while(<INFILE>){ 
  chomp;
  next if /^#/;
  my ($score,$name1,$strand1,$size1,$zstart1,$end1,$name2,$strand2,$size2,$zstart2,$end2,$nmatch,$nmismatch,$ngap,$continuity,$conPct,$identity,$idPct,$coverage,$covPct)=split(/\t/);
  my $algnBases=$coverage;
  $algnBases=~s/\/[0-9]+//;
  $lengths{$name1}=$size1;
  $lengths{$name2}=$size2;
  $largestContig=$name1;
  $smallestContig=$name2;
  if ($size2 > $size1){ die "The second contigs should be the smaller one.\n";}
  $idPct=~s/%//;
  $covPct=~s/%//;
  if($covPct >= $minHSPcov && $idPct >= $minHSPident){
#  if(($algnBases >= $minHSPAlgnBases || $covPct >= $minHSPcov) && $idPct >= $minHSPident){
#  if($algnBases >= $minHSPAlgnBases && $idPct >= $minHSPident){
   push @ranges, Bio::Range->new(-start=>$zstart2, -end=>$end2, -strand=>$strand2);
  }
 }
 my $nonOverlapingRanges='';
 $nonOverlapingRanges=reduceRanges(\@ranges,1,\%rangeSeen);
 my $totalLenAlgn=0;
 foreach my $r(@$nonOverlapingRanges){
#  print "$smallestContig ".$r->toString()."\n";
  my $lenAlgn= $r->end - $r->start + 1;
  $totalLenAlgn=$totalLenAlgn+($lenAlgn);
 }
 my $fractionLenAlgn=sprintf("%.2f", ($totalLenAlgn*100)/$lengths{$smallestContig});
 if ($fractionLenAlgn >= $minTotalFracAlgn){
  print CC4R "$largestContig\t$lengths{$largestContig}\t$smallestContig\t$lengths{$smallestContig}\t$totalLenAlgn\t$fractionLenAlgn\n";
  print "$largestContig\t$lengths{$largestContig}\t$smallestContig\t$lengths{$smallestContig}\t$totalLenAlgn\t$fractionLenAlgn\n";
  $contigSet1{$largestContig}=1;
  $contigSet2{$smallestContig}=1;
 }
}
close CC4R;

foreach my $contig(keys %contigSet2){
 if($contigSet1{$contig}){
  print STDERR "Contig already in largest contigs list, removing from candites for removal.\n";
 }
 else{
  print C2H "$contig\n";
 }
}
close C2H;

sub reduceRanges{
 my $ran=shift;
 my @ran=@$ran;
 my $overlapCount=shift;
 my @ranNew=();
 my $overlapsPresent=0;
 my $rangeSeen=shift;
 my %rangeSeen=%$rangeSeen;
 if($overlapCount>0){
  for (my $i=0;$i<@ran;$i++){
   my $overlapCountNew=0;
   next if $rangeSeen{$ran[$i]->start}{$ran[$i]->end};
   for (my $j=$i+1;$j<@ran;$j++){
    next if $rangeSeen{$ran[$j]->start}{$ran[$j]->end};
    if($ran[$i]->overlaps($ran[$j])){
     $overlapCountNew++;
     $overlapsPresent++;
     my ($start, $stop, $strand) = $ran[$i]->union($ran[$j]);
     push @ranNew, Bio::Range->new(-start=>$start, -end=>$stop, -strand=>$strand);
     $rangeSeen{$ran[$j]->start}{$ran[$j]->end}=1;
    }
   }
   if($overlapCountNew==0){
    push @ranNew, Bio::Range->new(-start=>$ran[$i]->start, -end=>$ran[$i]->end, -strand=>$ran[$i]->strand);
   }
  }
  reduceRanges(\@ranNew,$overlapsPresent,\%rangeSeen);
 }
 else{
  return \@ran
 }
}
