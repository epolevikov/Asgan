# Description
Asgan – [As]sembly [G]raphs [An]alyzer – is a tool for analysis of assembly graphs.
The tool takes two assembly graphs in the GFA format as input and finds the minimum set of
homologous sequences (synteny paths) for the graphs and then calculates different statistics
based on the found paths.

# Installation
```
git clone --recurse-submodules https://github.com/epolevikov/Asgan
make -C Asgan/lib/minimap2
```

# Usage
To run Asgan for the graphs in the _test_ folder, use the following command:
```
cd Asgan
python asgan.py \
    --input-query=test/NCTC9016-Flye.gfa \
    --input-target=test/NCTC9016-Canu.gfa \
    --out-dir=NCTC9016-Flye-vs-Canu
```
After analysis is finished, an output folder will contain:
* adjacency_graph_{query, target}.gv – visualization of synteny paths for the graphs.
* _synteny_paths.txt_ – synteny paths for the graphs.
* _stats.txt_ – different statistics based on the found synteny paths.

For the graphs above, synteny paths look like this:

<p align="center">
    <img src="https://github.com/epolevikov/Asgan/blob/master/graph-examples/flye_vs_canu.png">
</p>

The graph built by Canu consists of one connected component. Two paths represent forward (-2, -1, -4, +3) and
reverse complement (-3, +4, +1, +2) strands of a bacterial chromosome. The graph built by Flye also consists of
one connected component, but two complementary paths are merged through common unresolved repeats. Although the
structures of the graphs are different, they share one synteny path that corresponds to a bacterial chromosome.

# WABI Supplementary

https://zenodo.org/record/3198701
