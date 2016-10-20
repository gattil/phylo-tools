#!/usr/bin/env python3

import sys
import findtrunk
import unittest
import dendropy

class Data():
    DATA = "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);"
    # Get tree input
    TREE_NEWICK = dendropy.Tree.get(data=DATA,
                             schema="newick")


class TestEdges(unittest.TestCase):

    def test_EdgeCount_NotZero(self):
        tmp_cedge_inst = findtrunk.CEdges()
        tmp_cedge_inst.count_edges(Data.TREE_NEWICK)
        self.assertGreater(tmp_cedge_inst.edge_number, 0)

    def test_EdgeCount_ClassInt(self):
        tmp_cedge_inst = findtrunk.CEdges()
        tmp_cedge_inst.count_edges(Data.TREE_NEWICK)
        self.assertIsInstance(tmp_cedge_inst.edge_number,int)


# Main execution routine
if __name__ == "__main__":
    unittest.main()
