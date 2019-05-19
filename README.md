# Project description
Asgan (ASsembly Graphs ANalyzer) is a tool for analysis of assembly graphs.
Asgan takes two assembly graphs in GFA format as input and finds the minimum set of
homologous sequences (synteny paths) for the graphs and calculates different
statistics based on the found paths.

# Installation
```
git clone https://github.com/epolevikov/Asgan
git clone https://github.com/lh3/minimap2 Asgan/lib/minimap2
make -C Asgan/lib/minimap2
```

# Usage
```
python -m asgan \
    --input-query=test/NCTC9016/graph_flye.gfa \
    --input-target=test/NCTC9016/graph_canu.gfa \
    --out-dir=NCTC9016-Flye-vs-Canu
```

A file named synteny_paths.txt contains synteny paths for the graphs,
a file stats.txt contains different statistics based on the found paths.
Two file in .gv format contain a visualization for the paths.
