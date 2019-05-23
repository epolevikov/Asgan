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

| <img src="https://github.com/epolevikov/Asgan/blob/master/graph-examples/flye.svg" width=250> | <img src="https://github.com/epolevikov/Asgan/blob/master/graph-examples/canu.svg" width=175> |
| *Flye* | *Canu* |

# WABI Supplementary

https://zenodo.org/record/3198701
