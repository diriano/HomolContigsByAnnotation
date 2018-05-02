Find Homologous Contigs By Gene Annotation
===========

python script to identify genes that are present on multiple contigs


Introduction
============

Genome assemblies of non-inbred diploid individuals may result in the assembly of "homologous contigs", contig pairs that capture homologous genomic regions but are sufficiently divergent from each other as to be assembled independently. This phenomenom often results in a genome length that is longer than expected. Unlike other tools which merge homologous regions based on DNA sequence similarity ([Redundans](http://nar.oxfordjournals.org/content/44/12/e113), [HaploMerger](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3409271/pdf/1581.pdf)), my method uses genome annotation information to identify pairs of homologous contigs. The source of the genome annotation is flexible but I recommend a conservative approach based on identifying conserved core (single-copy) genes using [BUSCO](http://busco.ezlab.org/). Contigs that lack BUSCO genes will be missed by this analysis so adding an transcript-based annotation is more sensitive, but also prone to error due to true gene duplications or ambiguous mapping of transcripts.

This approach has been effective at deduplicating primary contigs in the Aedes aegypti FALCON-Unzip assembly, see my [poster](http://www.pacb.com/wp-content/uploads/Kingan-AGBT-2017-A-high-quality-genome-assembly-of-SMRT-sequences-diploid-mosquito-Aedes-aegypti.pdf) presented at AGBT 2017.

User Inputs
===========

At the command line, the user must specify the file paths for the following inputs:

1. BED file containing annotated genes and their gene spans (do not include mRNA, exons etc!). The shell script [BUSCOfull2BED](https://github.com/skingan/HomolContigsByAnnotation/blob/master/BUSCOfull2BED.sh) will produce this for duplicated genes.
2. Tab delimited text file with contig IDs and lengths.

Output
===========

Tab delimited text file with the following fields:

1. contigA ID
2. contigB ID
3. contigA length (always the longer contig) 
4. contigB length
5. number of genes shared for contig pair
6. comma-delimited list of shared gene IDs 


Downstream analysis recommendations
===========

A non-redundant list of contigB IDs is a reasonable starting list of contigs that should be flagged as allelic (haplotype) variants. Simply substracting their total length from the genome assembly length may result in a more accurate estimate of genome size.

_NOTE: This approach has been fruitful for PacBio assemblies with contig N50 stats on the order of > 1Mbp. I cannot vouch for its effectiveness on more fragmented NGS-bases assemblies._

Extension (DMRP)
===========

The original HomolContigsByAnnotation, reported a list of contig pairs that contain duplicated BUSCO. My extension carries out a LASTZ alignment between the two contigs and plots a dotplot (using R's ggplot2), highligting the position of identified BUSCO genes. This extension includes two scripts: compareContigPair.pl (requires BioPerl) and compareContigPairPlotAlignment.Rscript. 

Usage:

compareContigPair.pl HomolContigsByAnnotation_res.txt contigs.fixLineLength.fasta full_table_curated.BUSCO.tsv 


the contigs  file (contigs.fasta) must be well formated. Just to make sure you can use EMBOSS's seqret:

seqret -auto contigs.fasta contigs.fixLineLength.fasta

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.








