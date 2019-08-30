Asgan – [As]sembly [G]raphs [An]alyzer – is a tool for analysis of assembly graphs.
The tool takes two assembly graphs in the GFA format as input and finds the minimum set
of homologous sequences (synteny paths) shared between the graphs. As output, Asgan
produces various statistics and a visualization of the found paths in the .gv format.

# Installation
```
git clone --recurse-submodules https://github.com/epolevikov/Asgan
make -C Asgan/lib/minimap2
```

# Usage example
The _test_ folder contains two bacterial assembly from the NCTC collections produced by Flye
and Canu assemblers. To run Asgan for these assemblies, use the following command:
```
cd Asgan
python asgan.py \
    --input-query=test/flye-nctc9016.gfa \
    --input-target=test/canu-nctc9016.gfa \
    --out-dir=flye-vs-canu
```
After analysis is finished, the output directory will contain:
* adjacency_graph_{query, target}.gv – a visualization of synteny paths for the graphs.
* _synteny_paths.txt_ – synteny paths in the .txt format.
* _stats.txt_ – various statistics based on the found synteny paths.

Here is how the visualization looks like:

<p align="center">
    <img src="https://github.com/epolevikov/Asgan/blob/master/example.png">
</p>

The graph built by Canu consists of one connected component. Two paths represent forward (+1, +2, +3, +4) and
reverse complement (-4, -3, -2, -1) strands of a bacterial chromosome. The graph built by Flye also consists of
one connected component, but two complementary paths are merged through common unresolved repeats. Although the
structures of the graphs are different, they share one synteny path that corresponds to a bacterial chromosome.

Description of _synteny_paths.txt_ [TBA]
Description of _stats.txt_ [TBA]

# Tuning alignment parameters

To find an alignment between two assemblies, _Asgan_ utilizes _minimap2_. By default, _minimap2_ is
used with the _asm10_ preset. If two species diverge much (sequence divergence >10%), _minimap2_
will not find any alignments between them. If your datasets diverge much, use either _map-pb_ or _map-ont_ preset
using _--minimap-preset_ argument.

# WABI Supplementary

https://zenodo.org/record/3198701
