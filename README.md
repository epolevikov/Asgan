# Description
Asgan -- [As]sembly [G]raphs [An]alyzer -- is a tool for analysis of assembly graphs.
The tool takes two assembly graphs in the GFA format as input and finds the minimum set of
homologous sequences (synteny paths) for the graphs and then calculates different statistics
based on the found paths.

# Installation
```
git clone --recurse-submodules https://github.com/epolevikov/Asgan
make -C Asgan/lib/minimap2
```

# Usage
```
python asgan.py \
    --input-query=test/NCTC9016-Flye.gfa \
    --input-target=test/NCTC9016-Canu.gfa \
    --out-dir=NCTC9016-Flye-vs-Canu
```

A file named _synteny_paths.txt_ contains synteny paths for the graphs,
a file _stats.txt_ contains different statistics based on the found paths.
Two file with _.gv_ format contain a visualization for the paths.

# WABI Supplementary

https://zenodo.org/record/3198701
