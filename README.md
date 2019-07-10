# Cluster Specific Split-read Finder


Algorithm for identification of reads overlaping breakpoint junctions - Split-reads.
This tool was specificly design for long insert genome sequencing libraries (liGS), with small reads (<50bp).

## Before starting

To be able to run this tool, the user must have:

+ Aligned the liGS files against the reference genome using [BWA alignment software](http://bio-bwa.sourceforge.net/), obtaining a **SAM alignment file**
+ Identified the **breakpoint region**, using one of the improper-pair clustering methods available (i.e. SVdetect, readPairCluster...):

The user must have diferent sets of coordinates, depending on the type of alteration to be analyzed:

![alt text](https://cld.pt/dl/download/0997c4b6-d3a2-4a92-b8ad-81063eac74aa/Esquema%20split%20reads%20%281%29.jpg "Types of alterations and inputs examples")

+ large green and yellow arrows depict the locations of the clusters
+  A and B indicate the chromosomes
+ a1, a2, a3, a4, b1, b2, b3, b4 indicate the regions where the improper read pairs that identify the breakpoint region are mapped

Using these coordinates, the algorithm will try and map the split reads, assuming that:

+ split-reads were not mapped in previous steps;
+ the pair of the split-read is mapped within the translocation cluster and is marked in the SAM file as “unmapped-mate”.


The algorithm can be divided in two parts: **data selection**, where the potential split-reads and breakpoint regions are selected and prepared for analysis; and **read processing**, where BWA tries to map iteratively the potentially split-reads against the breakpoint regions.

**Data selection consists of:**

1. retrieving the FASTA sequences of the narrowest breakpoint intervals through NCBI API and BWA indexing for posterior use as mapping reference;
2. selection, from the SAM file, of unmapped mate read-pairs localized within the above defined genomic regions, for posterior processing.

**The read processing includes:**

3. alignment of the first and last 5 bp, designated as read chunks, of each unmapped read, against the reference sequences, and storage of the mapping data;
4. repeated realignment of chunks after sequentially increasing their size by 1 bp until no read chunk has a possible alignment or until the chunk size reaches the unmapped mate read size;
5. validation of alignment results, outputting only those where read chunks of an unmapped mate were mapped to different breakpoint regions, and the sum of the length of the chunks is equal or greater to the unmapped-mate read.


## Dependencies:
+ python2
+ python sys, os, collections and [biopython](https://github.com/biopython/biopython)
+ [BWA alignment software](http://bio-bwa.sourceforge.net/) , installed and on path


## Usage:


In the command line:
<pre><code> python split_reads_V5_beta.py [A:a-a’-b’-b] [B:c-c’-d’-d] [SAM file] [reference genome] [trans/inv/del/dup/ins]
</code></pre>

Where:
+ **[chr:a-a’-b’-b]** are the cluster coordinates from chrA 
+ **[chr:c-c’-d’-d]** are the cluster coordinates from chrB
**The coordinates may vary depending on the type of alteration. see the figura above**
+ **[SAM file]** is the SAM alignment file
+ **[reference genome]** is hg19/hg38
+ **[trans/inv/del/dup/ins]** is the type of ateration to be analyzed: trans for translocation, inv for inversion, del for deletion, dup for duplication, ins for intra and interchromosomal insertion

## The output:

![alt text](https://cld.pt/dl/download/a0ba07da-1053-4120-8c81-befbfb7f92ac/resume.png "Split results example")

The results are outputed to the terminal prompt as the example above.
The result is a set of tab-separated columns, with the information:
+ **Name of the original read** - name of the read before being split into chucks
+ **Chunk read name** - the name of the read after being split into chuncks. its the same as "Name of the original read" with the sufix "_start", "_end", to indicate the part of the chunck
+ **Orientation** - reverse or forward, depending on which orientation the read was mapped
+ **Size** - size of the mapped chunk
+ **Position** - mapping position of the chunk
+ **Sequence** - Fasta sequence of the mapped chunk
+ **Mate position** - mapping position of the read pair that was correctly mapped previously (for comparison)
+ **Mate orientation** - mapping orientation of the read pair that was correctly mapped previously (for comparison)


## Notes:

+ The tool was tested for translocations, inversions, deletions, duplications and insertions. For complex rearrangements one may decompose the rearrangement into simpler ones.
+ The tool was developed for small reads, but can, be used for bigger reads, including mate-pair variable size reads
+ Depending on the size of the regions and the coverage of the library, this tool might take a while to run
+ The results are outputed to the terminal prompt. The user may change that by adding to the end of the comand "> split_read.results"
+ Results may be validated manually based on the plausibility of their genomic positions, their orientation relative to each other, to the karyotype and to the reciprocal breakpoint.


## License:

GPLv2


## Found a bug?

Or maybe just wanto to drop some feedback? Just open an issue on github!