# bmotif: counting motifs in bipartite networks

## Overview

`bmotif` is software to count occurrences of motifs in bipartite networks, as well as the number of times each node appears in each unique position within motifs. `bmotif` considers all 44 unique bipartite motifs up to six nodes and all 148 unique positions within these motifs. It is available in [R](https://github.com/SimmonsBI/bmotif), [MATLAB](https://github.com/SimmonsBI/bmotif-matlab) and [Python](https://github.com/SimmonsBI/bmotif-python). `bmotif` was originally developed to analyse bipartite species interaction networks in ecology but its methods are general and can be applied to any bipartite graph. The code is released under the MIT license.

## Motif and motif position dictionary
The motifs corresponding to each motif ID and the positions corresponding to each motif position ID can be found in **Simmons, B.I, Sweering, M. J. M, Dicks, L. V., Sutherland, W. J., Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. (full citation will be placed here)**. These were defined following Baker et al (2015) Appendix 1 Figure A27.

## Dependencies

`bmotif` is written in Python 2.7 and uses the following modules:
- numpy
- pandas

## Scripts
- **BMotifs.py**: contains functions for counting the occurrences of motifs and the number of times each node appears in each unique position within motifs
- **guide.ipynb**: Shows how to use the two main functions on a simple example network

## License
The code is released under the MIT license (see LICENSE file).

## Citation
If you use the package in your work, please cite:
Simmons, B.I, Sweering, M. J. M, Dicks, L. V., Sutherland, W. J., Di Clemente, R. bmotif: a package for counting motifs in bipartite networks. (full citation will be placed here)

## References
Baker, N.J., Kaartinen, R., Roslin, T. and Stouffer, D.B., 2015. Species’ roles in food webs show fidelity across a highly variable oak forest. Ecography, 38(2), pp.130-139.
