# bmotif: counting motifs in bipartite networks

## Overview

`bmotif` is software to count occurrences of motifs in bipartite networks, as well as the number of times each node appears in each unique position within motifs. `bmotif` considers all 44 unique bipartite motifs up to six nodes, within which there are 148 unique positions. It is available in [R](https://github.com/SimmonsBI/bmotif), [MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and [Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was originally developed to analyse bipartite species interaction networks in ecology but its methods are general and can be applied to any bipartite graph. The code is released under the MIT license.

## Dependencies

`bmotif` is written in Python 2.7 and uses the following modules:
- numpy
- pandas

## Scripts
- **BMotifs.py**: contains functions for counting the occurrences of motifs and the number of times each node appears in each unique position within motifs
- **guide.ipynb**: Shows how to use the two main functions on a simple example network
