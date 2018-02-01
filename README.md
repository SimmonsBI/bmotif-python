# bmotif: counting motifs in bipartite networks

## Overview

`bmotif` is software to count occurrences of motifs in bipartite networks, as well as the number of times each node appears in each unique position within motifs. It is available in [R](https://github.com/SimmonsBI/bmotif), [MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and [Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was originally developed to analyse bipartite species interaction networks in ecology but its methods are general and can be applied to any bipartite graph. The code is released under the MIT license.

## Use

`bmotif` considers all 44 unique bipartite motifs up to six nodes. Within these motifs there are 148 unique positions.

`bmotif` has two functions: `count_motifs` and `count_pos_motifs`. `count_motifs` counts occurrences of motifs in a bipartite network. `count_pos_motifs` counts the number of times each node in a network occurs in each of the positions within the motifs.

## Dependencies

`bmotif` is written in Python 2.7 and uses the following modules:
- numpy
- pandas

## Scripts
- **BMotifs.py**: contains functions for counting the occurrences of motifs and the number of times each node appears in each unique position within motifs
- **guide.ipynb**: Shows how to use the two main functions on a simple example network
