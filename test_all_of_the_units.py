#!/usr/bin/env python3

import sys
import findtrunk
import unittest
import dendropy

class TestEdges(unittest.TestCase):
    def test_EdgeCount(self):
        # Get tree input
        tree = dendropy.Tree.get(data="((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);",
                                 schema="newick")
        tmp_cedge_inst = findtrunk.CEdges()
        tmp_cedge_inst.count_edges(tree)
        self.assertGreater(tmp_cedge_inst.edge_number, 0)


# Main execution routine
if __name__ == "__main__":
    unittest.main()
