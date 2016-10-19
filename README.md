# phylo-tools [![Build Status](https://travis-ci.org/gattil/phylo-tools.svg?branch=master)](https://travis-ci.org/gattil/phylo-tools)
This repository contains tools for phylogenetic tree manipulation

## `findtrunk.py`

This script reconstructs the phylogenetic tree trunk using a reverse-traversal method. The tree
topology is traversed backward starting iteratively from each leaf in the tree. The program
counts how many times each edge is traversed and it stores the sum of the counts in a new label
called "trunk" associated to each edge in the topology.

        $ python3 findtrunk.py -i input_tree.tree --input-format nexus

[Goto documentation](http://lorenzogatti.me/phylo-tools/findtrunk.m.html)


## `trunktraitevolution.py`

This script computes the permanence and the number of switches of a discrete trait on  the tree
trunk.

        $ python3 trunktraitevolution.py -i input_tree.tree --input-format nexus --feature location

[Goto documentation](http://lorenzogatti.me/phylo-tools/trunktraitevolution.m.html)
