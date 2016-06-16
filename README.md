# phylo-tools
This repository contains tools for phylogenetic tree manipulation

## `findtrunk.py`

This script reconstructs the phylogenetic tree trunk using a reverse-traversal method. The tree
topology is traversed backward starting iteratively from each leaf in the tree. The program
counts how many times each edge is traversed and it stores the sum of the counts in a new label
called "trunk" associated to each edge in the topology.

        $ python3 findtrunk.py -i input_tree.tree --input-format nexus
        
        