#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""

This script reconstructs the phylogenetic tree trunk using a reverse-traversal method. The tree
topology is traversed backward starting iteratively from each leaf in the tree. The program
counts how many times each edge is traversed and it stores the sum of the counts in a new label
called "trunk" associated to each edge in the topology.

## Example


    $ python3 findtrunk.py -i input_tree.tree --input-format nexus


## Arguments:

#### mandatory arguments:

1. -i INPUT_FILE  Input filepath of the tree file
2. --input-format INPUT_FORMAT String indicating tree file format (i.e. nexus, newick, phy, ...)

#### optional arguments:

1.  -o/--out OUTPUT_FILE  Output tree file name
2.  -l/--label DATA_LABEL Additional label for the data contained in the tree file
3.  --log LOG_LEVEL       Log level for the routine (--log=INFO)
4.  --log_to_file LOG_TO_FILE Log filename for the routine

## Compatibility
`findtrunk` has been tested on Python 3.4


## Contributing
`findtrunk` [is on GitHub](https://github.com/gattil/phylo-tools). Pull requests and bug reports
are welcome.

## Licence
`findtrunk` is in the public domain under MIT licence

> The MIT License (MIT)

> Copyright (c) 2016 Lorenzo Gatti

> Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

> The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""
# Import section (built-in modules|third-party modules)
import sys
import argparse
import dendropy
import os.path
import logging

# Authorship information

__project__ = 'phylo-tools'
__product__ = 'findtrunk'
__editor__ = 'PyCharm'
__author__ = 'lorenzogatti'
__copyright__ = "Copyright 2016, Applied Computational Genomics Team (ACGT)"
__credits__ = ["Lorenzo Gatti"]
__license__ = "GPL"
__date__ = '16/06/16'
__version__ = "1.0"
__maintainer__ = "Lorenzo Gatti"
__email__ = "lorenzo.gatti@zhaw.ch"
__status__ = "Development"

# Classes

class CEdges:
    """
       This class contains all the methods associated to object of class Edge

    """
    def __init__(self):
        self.edge_number = int()

    def count_edges(self, tree):
        """
           This method counts all the edges in the tree topology

           Returns:
               It assign the total number of edges to the edge_number object in the self contenitor
        """
        self.edge_number = 0

        for edge in tree.levelorder_edge_iter():
            self.edge_number += 1




# Routines
def arg_parser():
    """
    This function parses the arguments passed from the command line to the script

    Returns:
        It returns an object containing all the recognised arguments properly formatted
    """

    parser = argparse.ArgumentParser(prog='findtrunk.py',
                                     description='Retrieve the phylogenetic trunk from the tree '
                                                 'topology using a reverse-tree-traversal '
                                                 'approach.')
    parser.add_argument("-i", "--in", type=str, dest="input_file", required=True,
                        help='Input tree file')
    parser.add_argument("--input-format", type=str, dest="input_format", required=True,
                        help='Input tree file format')
    parser.add_argument("-o", "--out", dest="output_file", type=str, required=False,
                        default='', help='Output tree file')
    parser.add_argument("-l", "--label", dest="data_label", type=str, required=False,
                        default='', help='Label for the data contained in the tree')
    parser.add_argument("--log", dest="log_level", type=str, required=False,
                        default='INFO', help='Log level for the routine (--log=INFO)')
    parser.add_argument("--log_to_file", dest="log_to_file", type=bool, required=False,
                        default='', help='Log filename for the routine')

    return parser.parse_args()


def main(args):
    """
    This function executes the routines required by the program to identify the tree trunk and it
    returns a fully labeled tree topology in a new tree file.

    `args` is a the object returned by `arg_parser` function.
    """

    # Prepare output variables
    if not args.output_file:
        args.output_file = os.path.dirname(args.input_file) + '/' + args.data_label + '_out.tree'

    # Prepare logging device
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log_level)

    logging.getLogger(__product__)
    if args.log_to_file:

        logging.basicConfig(filename=args.log_file,
                            level=numeric_level,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:

        args.log_file = os.path.dirname(args.input_file) + '/' + args.data_label + '_out.log'

        logging.basicConfig(level=numeric_level,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.debug(__project__+":"+__product__+" - Execution started")

    # Read the tree file
    tree = dendropy.Tree.get(path=args.input_file,
                             schema=args.input_format,
                             extract_comment_metadata=True)

    # Initialise feature-storing-objects
    edgedict = dict()
    edgelist = list()

    # ------------------------------------------------------------------------------------------
    # Add UIDs to each node/edge in the tree
    logging.debug(__project__ + ":" + __product__ + " - Adding UIDs to nodes")
    i = 0
    for node in tree.levelorder_node_iter():
        node.label = i
        sys.stdout.write("\r%s | Labelling node: %3d" % (os.path.basename(args.input_file), i))
        sys.stdout.flush()
        i += 1
    sys.stdout.write("\n")

    logging.debug(__project__ + ":" + __product__ + " - Adding UIDs to edges")
    i = 0
    for edge in tree.levelorder_edge_iter():
        edge.label = i
        edgedict[i] = 0
        sys.stdout.write("\r%s | Labelling edge: %3d" % (os.path.basename(args.input_file), i))
        sys.stdout.flush()
        i += 1
    sys.stdout.write("\n")

    # ------------------------------------------------------------------------------------------
    # Infer the trunk from tree topology using a reverse traversal (from leaf-to-root)
    logging.debug(__project__ + ":" + __product__ + " - Computing relative root distances")
    tree.max_distance_from_root()
    for leaf in tree.leaf_node_iter():

        backwardnode = leaf.parent_node

        while backwardnode.root_distance > 0:
            edgelist.append((backwardnode.edge.label, 1))

            backwardnode = backwardnode.parent_node

    for pair in edgelist:
        edgedict[pair[0]] += pair[1]

    # ------------------------------------------------------------------------------------------
    # Add feature on topology
    for edge in tree.levelorder_edge_iter():
        edge.annotations.add_new(name="trunk", value=edgedict[edge.label])

    # Count all the edges in the tree topology (TEST)
    inst_edges = CEdges()
    inst_edges.count_edges(tree)
    logging.debug(__project__ + ":" + __product__ + " - Total number of edges is " + str(inst_edges.edge_number))

    # ------------------------------------------------------------------------------------------
    # Save new tree topology
    tree.write(path=args.output_file, schema='nexus')
    logging.debug(__project__+":"+__product__+" - Featured tree saved in: " + args.output_file)


# Main execution routine
if __name__ == "__main__":

    # Parse calling arguments
    args = arg_parser()
    # Call main routine
    main(args)

    # Exit program
    sys.exit(0)
