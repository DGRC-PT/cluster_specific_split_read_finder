# Cluster Specific Split-read Finder


Algorithm for identification of reads overlaping translocation breakpoint junctions - Split-reads.
This tool was specificly design for long insert genome sequencing libraries (liGS), with small reads (<50bp).

## Before starting

To be able to run this tool, the user must have:

+ Aligned the liGS files against the reference genome using [BWA alignment software](http://bio-bwa.sourceforge.net/), obtaining a **SAM alignment file**
+ Identified the **breakpoint region**, using one of the improper-pair clustering methods available (i.e. SVdetect, readPairCluster...):

The user must have two sets of coordinates for each rearrangement breakpoint:

![alt text](https://cld.pt/dl/download/9b25e73b-b2ed-47ca-8dea-644442f6600d/clust.png "Read cluster example")

+  **a-a', b-b', c-c', d-d'** - indicate the regions where the improper read pairs that identify the breakpoint region are mapped (translocation cluster)
+ **a'-b', c'-d'** - indicate the breakpoint region

Using this coordinates, the algorithm will try and map the split reads, assuming that:

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
<pre><code> python split_reads_V5_beta.py [chr:a-a’-b’-b] [chr:c-c’-d’-d] [SAM file] [reference genome]
</code></pre>

Where:
+ **[chr:a-a’-b’-b]** are the cluster coordinates from chrA
+ **[chr:c-c’-d’-d]** are the cluster coordinates from chrB
+ **[SAM file]** is the SAM alignment file
+ **[reference genome]** is hg19/hg38


## The output:

![alt text](https://cld.pt/dl/download/84b984ac-e030-46f7-8c30-1afa9acef55b/splits.png "Split results example")

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

+ The tool was only tested for translocations, but may work with other type of SVs
+ The tool was developed for small reads, but can, theoreticly, be used for bigger reads
+ Depending on the size of the regions and the coverage of the library, this tool might take a while to run
+ The results are outputed to the terminal prompt. The user may change that by adding to the end of the comand "> split_read.results"
+ Results may be validated manually based on the plausibility of their genomic positions, their orientation relative to each other, to the karyotype and to the reciprocal breakpoint.


## License:

GPLv2


## Found a bug?

Or maybe just wanto to drop some feedback? Just open an issue on github!
   