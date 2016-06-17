#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""

This script computes the permanence and the number of switches of a discrete trait on  the tree
trunk.

## Example


    $ python3 trunktraitevolution.py -i input_tree.tree --input-format nexus --feature location


## Arguments:

#### mandatory arguments:

1. -i INPUT_FILE  Input filepath of the tree file
2. --input-format INPUT_FORMAT String indicating tree file format (i.e. nexus, newick, phy, ...)

#### optional arguments:

1.  -o/--out OUTPUT_FILE  Output tree file name
2.  -l/--label DATA_LABEL Additional label for the data contained in the tree file
3.  --feature FEATURE_ANNOTATION Discrete trait to monitor for trunk switches
4.  --trunk-threshold TRUNK_THRESHOLD Integer value indicating the lower bound to define the tree
                                      trunk. When not defined, all the tree will be considered
4.  --log LOG_LEVEL       Log level for the routine (--log=INFO)
5.  --log_to_file LOG_TO_FILE Log filename for the routine

## Compatibility
`trunktraitevolution` has been tested on Python 3.4


## Contributing
`trunktraitevolution` [is on GitHub](https://github.com/gattil/phylo-tools). Pull requests and bug
reports are welcome.

## Licence
`trunktraitevolution` is in the public domain under MIT licence

> The MIT License (MIT)

> Copyright (c) 2016 Lorenzo Gatti

> Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

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
import csv

# Authorship information

__project__ = 'phylo-tools'
__product__ = 'trunktraitevolution'
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
    parser.add_argument("--feature", dest="feature_annotation", type=str, required=True,
                        default='', help='Discrete trait')
    parser.add_argument("--trunk-threshold", dest="trunk_threshold", type=int, required=False,
                        default=0, help='Trunk value threshold')
    parser.add_argument("-l", "--label", dest="data_label", type=str, required=False,
                        default='', help='Label for the data contained in the tree')
    parser.add_argument("--log", dest="log_level", type=str, required=False,
                        default='INFO', help='Log level for the routine (--log=INFO)')
    parser.add_argument("--log_to_file", dest="log_to_file", type=bool, required=False,
                        default='', help='Log filename for the routine')

    return parser.parse_args()


def main(args):
    """
    This function executes the routines required by the program to identify the number of switches
    and the permanence of the discrete trait on the tree trunk

    `args` is a the object returned by `arg_parser` function.
    """

    # Prepare output variables
    if not args.output_file:

        filename = os.path.dirname(args.input_file) + '/' + args.data_label + '_' \
                   + args.feature_annotation

        args.output_file = filename

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

    logging.debug(__project__ + ":" + __product__ + " - Execution started")

    # Read the tree file
    tree = dendropy.Tree.get(path=args.input_file,
                             schema=args.input_format,
                             extract_comment_metadata=True)

    feature_permanence = dict()

    with open(args.output_file+'_switches.csv', 'w') as csvfile:

        fieldnames = ['FROM-ID', 'TO-ID', 'F-AGE', 'T-AGE', 'DURATION', 'VFROM', 'VTO', 'C']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        sc = 0
        # per each node in the tree topology, loop over the child nodes and retrieve the length
        # of the branch grouped according to
        for node in tree.preorder_node_iter():

            # check if the node has been already visited before
            if 'visited' not in node.annotations.values_as_dict().keys():

                # get the requested annotation for the parent node
                parent_annotation = node.annotations.values_as_dict()[args.feature_annotation]
                parent_id = node.label
                parent_height = node.annotations.values_as_dict()['height']

                # Per each child node in the tree starting from the node selected in
                # preorder-traversing
                for child in node.child_nodes():

                    # Check if the child has been labelled with a number
                    if child.label:
                        # count the number of switches for the discrete trait occurring on the
                        # trunk of the tree
                        if int(child.annotations.values_as_dict()['trunk']) > args.trunk_threshold:

                            # Get annotation of the current child node
                            child_annotation = child.annotations.values_as_dict()[
                                args.feature_annotation]

                            # Compute the permanence
                            if parent_annotation not in feature_permanence:
                                feature_permanence[parent_annotation] = {}

                            if child_annotation not in feature_permanence[parent_annotation].keys():
                                feature_permanence[parent_annotation][child_annotation] = \
                                    child.edge.length
                            else:
                                feature_permanence[parent_annotation][child_annotation] += \
                                    child.edge.length

                            #    feature_permanence[child_annotation] += node.edge.length
                            # else:
                            #    feature_permanence[child_annotation] = node.edge.length

                            # Call switches

                            if parent_annotation != child_annotation:
                                c = 1
                                sc += 1
                            else:
                                c = 0

                            # Store
                            writer.writerow({'FROM-ID': parent_id,
                                             'TO-ID': child.label,
                                             'F-AGE': parent_height,
                                             'T-AGE': child.annotations.values_as_dict()['height'],
                                             'DURATION': child.edge.length,
                                             'VFROM': parent_annotation,
                                             'VTO': child_annotation,
                                             'C': c})

                            # Re-assigning internal values
                            parent_annotation = child_annotation
                            parent_id = child.label
                            parent_height = child.annotations.values_as_dict()['height']

                            # Complete visiting the node, adding an annotation indicating the
                            # successful visit
                            child.annotations.add_new(name="visited", value=1)

    logging.info(__project__ + ":" + __product__ + " - The discrete trait [" +
                 args.feature_annotation + '] shows ' + str(sc) + ' switches on the trunk')

    with open(args.output_file+'_summary.csv', 'w') as csvfile:

        fieldnames = ['VFROM', 'VTO', 'DURATION']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for vfrom in feature_permanence.keys():
            for vto in feature_permanence[vfrom].keys():
                if vfrom == vto:
                    writer.writerow({'VFROM': vfrom,
                                     'VTO': vto,
                                     'DURATION': feature_permanence[vfrom][vto]})

    # pprint.pprint(feature_permanence, width=1)

# Main execution routine
if __name__ == "__main__":

    # Parse calling arguments
    parsed_args = arg_parser()
    # Call main routine
    main(parsed_args)

    # Exit program
    sys.exit(0)
