![nanochrome](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/nanochrome3.svg)

Combines low depth 10x Chromium linkage with low depth nanopore long-reads to scaffold difficult assemblies

## Quickstart

No compilation required.

Clone the scaffolder:
```
git clone https://github.com/OliverCardiff/Nanochrome_scaffolder
```

Typical Usage:
```
./Run_Nanochrome.sh -g genome.fa -r chromium_reads.fq -n nanopore_reads.fq -f 29000 -p my_prefix -l 6
```

### Dependencies:


1. [Python3](https://realpython.com/installing-python/#ubuntu)

Follow the link and install python3, if you don't have it already.

2. [minimap2](https://github.com/lh3/minimap2)

Minimap2 can either downloaded as binaries on linux, or built from source. Follow the above link to find out how to do this.

Quickly grabbing the minimap2 binaries and added them to the PATH temporarily:
```
curl -L https://github.com/lh3/minimap2/releases/download/v2.16/minimap2-2.16_x64-linux.tar.bz2 | tar -jxvf -
cd minimap2-2.16_x64-linux/
export PATH=$PATH:$PWD
```

### Usage messages:

Shell script:
```
Usage: ./Run_Nanochrome.sh
        -g <file> Genome to scaffold in fasta format
        -r <file(,file)> Chromium paired reads (if one file will assume interleaved)
        -n <file> Nanopore reads, can be fasta, or fastq
        -f <integer> 10x Library mean fragment length
        -p <string> Prefix for output/processing files

        OPTIONAL:
        -l <integer> 3 - 20, Leniency: error/contiguity trade-off. Higher = More Error [5]
        -c <flag> (clean up) Add flag if you wish to clean up non-essential processing files
```

chrome_candidates.py:
```
Nanochrome - building barcode linkage candidates

ARG1: <genome.fa> A Genome file
ARG2: <alns.sam> Mapped 10X library
ARG3: <integer> Fragment Length (x10 library)
ARG4: <prefix> Unique prefix for outputs

Output:

<prefix>_candidates.fa: A set of scaffold candidates
<prefix>_table.tsv: Potential connection table
<prefix>_network.tsv: A network file for visualisation
<prefix>_nodes.tsv: A node description file for visualisation

You need to specify the three positional arguments
```

nano_confirms.py:
```
Nanochrome - confirming 10x with nanopore

Input:

ARG1: <genome.fa> A Genome file
ARG2: <nc_table.tsv> output from chrome_candidates.py
ARG3: <alns.paf> Mapped Longread library
ARG4: <integer> Fragment Length (x10 library)
ARG5: <integer> Tangle leniency [5-15]
ARG6: <prefix> Unique prefix for outputs

Output:

<prefix>_nc_scaffolded.fa: Final Scaffolded Genome
<prefix>_network_final.tsv: A network file for visualisation
<prefix>_nodes_final.tsv: A node description file for visualisation
```

## How It Works:

[Chromium linked reads](https://www.10xgenomics.com/technology/) from 10x genomics work by barcoding many reads from the same long molecules of DNA. However, the weakness of this technology is that the fragment -> barcode assignment is stochastic, and mutliple fragments of DNA are typically associated with a single barcode. With sufficient sequencing depth, the original molecule structures can still be resolved for most model organisms. However, these organisms are typically vertebrates, or super inbred creepy crawlies.

First, some form of de Brujin graph is used to assemble contigs, then the physical linkage information is used to build these into longer scaffolds, like so:

![scaff_idea](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/flow_chart.svg)

In the above image we can see that the multiple contigs aligned to by a single barcode's reads represent a highly ambiguous situation for scaffolding. However, the probability of very very many barcodes linking the same two scaffolds by chance is so incredibly low that we are safe to assume this means they are actually physically connected!


But what if you have an earthworm with 10% absolute base sequence divergence between alleles? Will it still work the same? What if your novel organism is absolutely riddled with phages and transposons!? In these cases the initial kmer graph will often not resolve long enough contigs to scaffold properly.. Why?

![doesntwork](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/can_scaff.svg)

The sequencing of fragments does not preserve their sense, and with the sparsity of coverage it is very difficult to orient by strand, or arrange by order, a series smaller fragments. Lets look at a different solution to the same problem...

[Oxford Nanopore](https://nanoporetech.com/learn-more), and specifically the minION sequencer, are capable of creating very long reads of single DNA molecules (up to, and even above 100kb). These reads can be very useful for scaffolding in cases where the mapping is non-redundant. In assemblies of low allelic divergence genomes, or those without extensive repeat content, simple repeated overlaps of paired contigs are sufficient to create well oriented scaffolds - and even close to accurate consensus sequence for the gaps between them.

![nanopores2](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/nanopore2.svg)

However, the error rate for each read is ~90-95% accuracy and this can create issues. It means that for each read, approximately 1-in-10 to 1-in-20 bases will be miscalled. If the absolute base sequence divergence of two alleles in an inflated assembly is comparable to this, it will introduce substantial ambiguity to the scaffolding process. For example:

![nanopores1](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/nanopore.svg)

From the above illustration it can be seen that there are assembly 'grey areas' in the construction of genomic scaffolds when using either of the two library types. By combining them in a single scaffolding process, Nanochrome seeks to allow for the the strengths of each to account for the other's flaws!

## The sort of thing chrome_candidates.py produces:
![demoimg](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/feature_img.png)
## The sort of thing nano_confirms.py produces:
![Close](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/test_zoomout.png)
## As close as you can get to the sort of thing nano_confirms.py produces without getting your eyes wet:
![Close](https://github.com/OliverCardiff/Nanochrome_scaffolder/blob/master/media/net_closeup.png)
