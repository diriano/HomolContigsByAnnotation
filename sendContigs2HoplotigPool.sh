#!/bin/bash

CONTIGS_FILE=../curated.fixLineLength.fasta
HAPLOTIG_FILE=../curated.haplotigs.fasta
HAPLOTIG_FILE_FINAL=haplotigs.final.fasta
CONTIGS2HAPLOTIGS_LIST=Contigs2HaplotigPool.txt
CONTIGS2PRIMARY_LIST=primary.final.txt
CONTIGS2PRIMARY=primary.final.fasta
#Fix Line length on haplotig file
seqret $HAPLOTIG_FILE $HAPLOTIG_FILE_FINAL

#Prepare list of contigs to be sent to the primary assembly file
grep ">" $CONTIGS_FILE | grep -v -f $CONTIGS2HAPLOTIGS_LIST  |cut -f 1 -d' ' |sed 's/>//' > $CONTIGS2PRIMARY_LIST
sed -i "s#^#$CONTIGS_FILE:#" $CONTIGS2PRIMARY_LIST
#Prepare list of contigs to be sent to the haplotig file
sed -i "s#^#$CONTIGS_FILE:#" $CONTIGS2HAPLOTIGS_LIST

#Extract contigs and send them to haplotig pool
seqret -auto @$CONTIGS2HAPLOTIGS_LIST stdout >> $HAPLOTIG_FILE_FINAL


#Extract contigs and send them to final primary assembly
seqret -auto @$CONTIGS2PRIMARY_LIST stdout > $CONTIGS2PRIMARY
