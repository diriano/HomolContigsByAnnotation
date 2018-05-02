#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::DB::Fasta;

die if !@ARGV;

my $pairsFile    =$ARGV[0];
my $sequenceFile =$ARGV[1];
my $genePosFile = $ARGV[2];
my $db= Bio::DB::Fasta->new($sequenceFile);
open PAIRS, $pairsFile;
while(<PAIRS>){
 chomp;
 next if /^Contig1\tContig2/;
 my ($contig1,$contig2,$len1,$len2,$sharedGenes,$geneIDs)=split(/\t/);
 my $outAlgn=$contig1."__".$contig2."_lastz.rdotplot";
 my $contig1File=$contig1.'.fasta';
 my $contig2File=$contig2.'.fasta';
 if(!-s $contig1File){
  my $contig1Seq     = $db->get_Seq_by_id($contig1);
  my $seq_out = Bio::SeqIO->new( -file   => ">$contig1File",
                                 -format => 'fasta',
                               );
  $seq_out->write_seq($contig1Seq);
 }
 if(!-s $contig2File){
  my $contig2Seq     = $db->get_Seq_by_id($contig2);
  my $seq_out = Bio::SeqIO->new( -file   => ">$contig2File",
                                 -format => 'fasta',
                               );
  $seq_out->write_seq($contig2Seq);
 }
 my $outGenesContig1=$contig1.".genePos.txt";
 my $outGenesContig2=$contig2.".genePos.txt";
 system("grep $contig1 $genePosFile > $outGenesContig1");
 system("grep $contig2 $genePosFile > $outGenesContig2");
 system("lastz $contig1File $contig2File --notransition --step=20 --gapped --format=rdotplot > $outAlgn");
 system("Rscript --vanilla /data1/bioinfo/HomolContigsByAnnotation/compareContigPairPlotAlignment.Rscript $outAlgn $outGenesContig1 $outGenesContig2");
 unlink $contig1File;
 unlink $contig2File;
 unlink $outAlgn;
 unlink $outGenesContig1;
 unlink $outGenesContig2;
}
close PAIRS;
