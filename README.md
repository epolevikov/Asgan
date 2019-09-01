Asgan – <strong>As</strong>sembly <strong>G</strong>raphs <strong>An</strong>alyzer – is a tool for analysis of
assembly graphs. The tool takes two assembly graphs in the _GFA_ format as input and finds the minimum set
of homologous sequences (synteny paths) shared between the graphs. As output, Asgan produces various statistics and
a visualization of the found paths in the _gv_ format.

# Installation
```
git clone --recurse-submodules https://github.com/epolevikov/Asgan
make -C Asgan/lib/minimap2
```

# Usage example
The _test_ folder contains two bacterial assembly from the NCTC collection produced by Flye
and Canu assemblers. To run Asgan for the datasets, use the command below:
```
cd Asgan
python asgan.py \
    --input-query=test/flye-nctc9016.gfa \
    --input-target=test/canu-nctc9016.gfa \
    --out-dir=flye-vs-canu
```
After analysis is finished, the output directory will contain the following files:
* <strong>adjacency_graph_{query, target}.gv</strong> – a visualization of synteny paths.
* <strong>synteny_paths.txt</strong> – synteny paths in the format of an alignment.
* <strong>stats.txt</strong> – various statistics for the graphs.

## Visualization

Here is how the visualization looks like:

<p align="center">
    <img src="https://github.com/epolevikov/Asgan/blob/master/example.png">
</p>

The graph built by Canu consists of two separated sequences. One of them represents a forward (+1, +2, +3, +4) strand,
the other corresponds to a reverse complement (-4, -3, -2, -1) strand of a bacterial chromosome. The graph built by Flye
consists of one connected component, where two complementary strands are merged through common unresolved repeats.
Although the structures of the graphs are different, they share one synteny path that corresponds to a bacterial chromosome.

## Synteny paths

A file named _synteny_paths.txt_ contains the found synteny paths in the format of an alignment. For the above
datasets, the file looks like this:
```
+1      contig_8+       5'078'954       56'910          389'537         contig_2-       425'024         47'878          378'220     
        contig_8+       5'078'954       389'537         451'517         contig_7+       16'794          0               16'794      
+2      contig_8+       5'078'954       451'517         662'367         contig_6+       209'014         0               209'014     
        contig_8+       5'078'954       662'367         682'556         contig_7-       16'794          0               16'794      
+3      contig_8+       5'078'954       682'556         3'667'386       contig_1+       2'964'866       11              2'964'858   
        contig_8+       5'078'954       3'667'386       3'691'634       contig_5-       24'049          0               24'049      
+4      contig_8+       5'078'954       3'691'634       5'064'579       contig_3+       1'364'661       0               1'364'661
```
The first column contains the names of the alignment blocks. The rows with ids (+1, +2, +3, +4) correspond to the
unique alignment blocks, while the rows with the empty block title represent repeats. The second column corresponds
to the sequences of the query assembly. The following three columns show the length of a sequence, the starting and
the ending position of an alignment block accordingly. The remaining columns correspond to the sequences, lengths,
and mapping positions of the alignment blocks for the target assembly.

## Statistics

Here is the content of a file named _stats.txt_:
```
        Query           Target
cc      1               1           
useqs   1               4           

seqs    1               6           
tlen    5'078'954       5'004'408   
N50     5'078'954       2'964'866   
L50     1               1           

blocks  4               4           
tlen    4'901'252       4'868'864   
bcvg    0.965           0.973       
N50     2'984'830       2'964'847   
L50     1               1           

paths   1               1           
tlen    5'007'669       4'926'501   
N50     5'007'669       4'926'501   
L50     1               1
```
The file consists of three columns: the first one contains the titles of various statistics, the remaining two show the
values of these statistics for the _Query_ and the _Target_ assemblies accordingly. The statistics are split into four
groups. The first group contains two statistics:

* __cc__ – the number of connected components in an assembly graph. Only the components containing at least one
alignment block are counted. Note that even if two complementary strands are separated from each other, they are
assumed to represent one component. Thus, although the graph built by Canu consists of two separated strands, the
number of connected components is equal to 1.
* __useqs__ – the number of unique sequences in an assembly graph. The number is calculated as the difference between the
total number of sequences and the number of repeats. Note that only the sequences that belong to a connected component
with at least one alignment block are counted. A sequence is considered to be a repeat either if its length is lower
than 50k or it does not contain alignment blocks.

The following three groups start with the statistics named __seqs__, __blocks__, and __paths__. Below is their
description:

* __seqs__ – the total number of sequences in an assembly graph. Only the sequences that belong to a connected component
with at least one alignment block are counted.
* __blocks__ – the total number of alignment blocks between two assembly graphs.
* __paths__ – the number of synteny paths shared between the assembly graphs. The paths are obtained
by merging the alignment blocks that appear in the same order both in the query and in the target graph. See the
manuscript [TBA] if you are interested in a more detailed description of the algorithm.

Important that for the last four statistics, two complementary sequences (contig_1+ and contig_1- for example)
are counted as one.

Each of the last three groups (__seqs__, __blocks__, __paths__) contains statistics named __tlen__ (total length),
__N50__, __L50__ that can be used to estimate the contiguity of the corresponding sequences. The __blocks__ group also
has a statistic named __bcvg__ (block coverage). This number is calculated as the total length of alignment blocks
divided by the total length of sequences.

# Use cases

## Assessing the quality of an assembly graph

When we run an assembler for a collection of reads, we expect the assembler to reconstruct each chromosome of a genome.
Thus, if a genome contains _N_ chromosomes, the ideal assembly graph would consist of _N_ connected components each
representing one chromosome. Due to repeats and errors in reads, assembly graphs often contains defects. For example,
the parts of different chromosomes might be merged in one connected component. Or, one chromosome might be separated
into several parts belonging to different connected components.

If you have a reference genome available, you may use Asgan to estimate how far or close your assembly graph from the
ideal one. To do that, run Asgan using the assembly graph as _Query_ and the reference as _Target_. The difference
between the number of chromosomes and the number of synteny paths will reveal the quality of the assembly graph.

Below you can see the results of such an analysis for the assembly graph built by Flye for the C.elegans dataset. First,
let's have a look at the statistics:
```
        Query           Target
cc      2               6           
useqs   78              6           

seqs    165             6           
tlen    98'905'858      100'272'607
N50     1'919'214       17'493'829  
L50     18              3           

blocks  80              80          
tlen    96'408'233      97'366'606  
bcvg    0.975           0.971       
N50     1'888'725       1'875'295   
L50     18              18          

paths   16              16          
tlen    97'713'780      99'613'824  
N50     8'650'665       8'851'958   
L50     4               4
```
The reference genome (_Target_) consists of six connected components (one for each chromosome), while the assembly
graph (_Query_) has only two. This tells us that the graph contains unresolved repeats and the parts of several
chromosomes are merged through them in one component. The number of common alignment blocks between the assembly graph
and the reference is 80, while the number of synteny paths is 16. This means, on the one hand, that most of the
alignment blocks appear in the assembly graph in the same order as in the reference indicating the correct structure of
the graph. On the other hand, sixteen is greater than 6 (the expected number of common paths), which indicates that the
graph contains some defects.

We don't provide a visualization for the assembly graph here since its structure is complicated and it is hard to draw
any conclusions just looking at it. Instead, it might be useful to see how the synteny paths traverse the reference
chromosomes:

<p align="center">
    <img src="https://github.com/epolevikov/Asgan/blob/master/chromosomes.png">
</p>

Synteny paths shared between the assembly graph a the reference are shown in different colors. You can see that one
chromosome is covered by one path (shown in blue), which is the ideal case expected for each of the chromosomes. The
remaining ones, however, are covered by two or more paths indicating that some parts of the graphs have defects. The
nodes incident to edges with different colors correspond to the defective parts of the assembly graph.

## Comparing two assembly graphs built by different assemblers

Assume that we applied Asgan for two fragmented bacterial assembly and got one synteny path for them. Thus,
having two fragmented assembly, we obtained a complete one using synteny paths decomposition. At this point,
one possible application of synteny paths is improving assembly quality. If two alignment blocks appear in the same
order in two different assembly graphs, these blocks are likely to appear in the same order in a reference genome.
In other words, a collection of synteny paths for two fragmented assembly graphs might represent more contiguous
sequences comparing with the initial sequences for each of the graphs.

We didn't performed thorough experiments based on this idea and it is likely that for real datasets one may encounter
some pitfalls. This, probably, might be a direction for the further research.

Link to the manuscript: TBA.

## Comparing assemblies of different species

Asgan can be used to compare assemblies of different species. We showed that the N50 metric for synteny paths correlates
with the genomic distance between species. See the manuscript for details.

Link to the manuscript: TBA.

# Tuning alignment parameters

To find an alignment between two assemblies, Asgan utilizes minimap2. By default, minimap2 is
used with the _asm10_ preset. In some cases, the preset might need to be changed. For example, if two assemblies
diverge much (sequence diverge > 10%), minimap2 will not find alignment blocks between them. For highly diverged
species, we recommend to use either _map-pb_ or _map_ont_ preset. The default preset can changed using the
_--minimap-preset_ argument.

# WABI Supplementary

https://zenodo.org/record/3198701
