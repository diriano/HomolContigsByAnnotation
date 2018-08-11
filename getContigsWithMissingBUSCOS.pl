#!/usr/bin/perl

use strict;
use warnings;
#runlike $0 /data2/Vellozia/CANU/Vintermedia/FirstTry/BUSCO/run_VintermediaFirstTryCANU.contigs.fasta.BUSCO_v3.0_plants/full_table_VintermediaFirstTryCANU.contigs.fasta.BUSCO_v3.0_plants.tsv /data2/Vellozia/CANU/Vintermedia/FirstTry/PurgeHaplotypes/BUSCO/run_curated.BUSCO_v3.0_plants/full_table_curated.BUSCO_v3.0_plants.tsv 
open BUSCO1, $ARGV[0];
open BUSCO2, $ARGV[1];
