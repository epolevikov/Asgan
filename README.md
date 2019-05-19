# Project description
Asgan (ASsembly Graphs ANalyzer) -- is a tool for analysis of assembly graphs.
Asgan takes two assembly graphs in GFA format as input and finds the minimum set of
homologous sequences for the graphs and calculates different statistics based
of the found paths.

# Installation
```
git clone https://github.com/epolevikov/Asgan
git clone https://github.com/lh3/minimap2 Asgan/lib/minimap2
make -C Asgan/lib/minimap2
```
