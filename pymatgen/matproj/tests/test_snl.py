#!/usr/bin/env python

"""
Unit tests for StructureNL (SNL) format
"""

from __future__ import division

__author__ = "Anubhav Jain"
__credits__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "2/14/13"

import datetime
import unittest
import numpy as np
import json

from pymatgen import Structure
from pymatgen.matproj.snl import StructureNL, HistoryNode
from pymatgen.serializers.json_coders import PMGJSONEncoder, PMGJSONDecoder


class StructureNLCase(unittest.TestCase):
    def setUp(self):
        # set up a Structure
        self.s = Structure(np.eye(3, 3) * 3, ["Fe"], [[0, 0, 0]])
        self.s2 = Structure(np.eye(3, 3) * 3, ["Al"], [[0, 0, 0]])

        # set up BibTeX strings
        self.matproj = "@misc{MaterialsProject,\ntitle = {{Materials " \
                       "Project}},\nurl = {http://www.materialsproject.org}\n}"
        self.pmg = "@article{Ong2013,\n author = {Ong, " \
                   "Shyue Ping and Richards, William Davidson and Jain, " \
                   "Anubhav and Hautier, Geoffroy and Kocher, " \
                   "Michael and Cholia, Shreyas and Gunter, Dan and Chevrier," \
                   " Vincent L. and Persson, Kristin A. and Ceder, Gerbrand}," \
                   "\n doi = {10.1016/j.commatsci.2012.10.028}," \
                   "\n issn = {09270256},\n journal = {Computational " \
                   "Materials Science},\n month = feb,\n pages = {314--319}," \
                   "\n publisher = {Elsevier B.V.}," \
                   "\n title = {{Python Materials Genomics (pymatgen): A " \
                   "robust, open-source python library for materials " \
                   "analysis}},\n url = {http://linkinghub.elsevier" \
                   ".com/retrieve/pii/S0927025612006295},\n volume = {68}," \
                   "\n year = {2013}\n}"
        repeat = "REPEAT" * 10000
        self.superlong = "@misc{SuperLong,\ntitle = {{" + repeat + "}}}"
        self.unicode_title = u"@misc{Unicode_Title,\ntitle = {{A \u73ab is a rose}}}".encode('utf_8')
        self.junk = "This is junk text, not a BibTeX reference"

        # set up some authors
        self.hulk = [{"name": "Hulk", "email": "hulk@avengers.com"}]
        self.america = "Captain America <captainamerica@avengers.com>"
        self.thor = [("Thor", "thor@avengers.com")]
        self.duo = "Iron Man <ironman@avengers.com>, " \
                   "Black Widow <blackwidow@avengers.com>"

        # set up HistoryNodes
        self.valid_node = HistoryNode("DB 1", "www.db1URLgoeshere.com",
                                      {"db1_id": 12424})
        self.valid_node2 = {"name": "DB 2", "url": "www.db2URLgoeshere.com",
                            "description": {"db2_id": 12424}}
        self.invalid_node = {"name": "DB 3",
                             "url": "http://www.db3isnotavalidnode.com"}

    def test_authors(self):
        a = StructureNL(self.s, self.hulk, references=self.pmg)
        self.assertEquals(a.authors[0].name, "Hulk")
        self.assertEquals(a.authors[0].email, "hulk@avengers.com")

        a = StructureNL(self.s, self.america, references=self.pmg)
        self.assertEquals(a.authors[0].name, "Captain America")
        self.assertEquals(a.authors[0].email, "captainamerica@avengers.com")

        a = StructureNL(self.s, self.thor, references=self.pmg)
        self.assertEquals(a.authors[0].name, "Thor")
        self.assertEquals(a.authors[0].email, "thor@avengers.com")

        a = StructureNL(self.s, self.duo, references=self.pmg)
        self.assertEquals(a.authors[0].name, "Iron Man")
        self.assertEquals(a.authors[0].email, "ironman@avengers.com")
        self.assertEquals(a.authors[1].name, "Black Widow")
        self.assertEquals(a.authors[1].email, "blackwidow@avengers.com")
        StructureNL(self.s, self.hulk, references=self.pmg)

    def test_references(self):
        # junk reference should not work
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          references=self.junk)

        # good references should be ok
        StructureNL(self.s, self.hulk, references=self.pmg)

        # unicode references should work
        StructureNL(self.s, self.hulk, references=self.unicode_title)

        # multi-line references should be OK
        StructureNL(self.s, self.hulk,
                    references='\n'.join([self.matproj, self.pmg]))

        # super long references are bad
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          references=self.superlong)

    def test_historynodes(self):
        a = StructureNL(self.s, self.hulk, history=[self.valid_node])
        self.assertEquals(a.history[0].name, "DB 1")
        self.assertEquals(a.history[0].url, "www.db1URLgoeshere.com")
        self.assertEquals(a.history[0].description, {"db1_id": 12424})

        a = StructureNL(self.s, self.hulk,
                        history=[self.valid_node, self.valid_node2])
        self.assertEquals(a.history[1].name, "DB 2")
        self.assertEquals(a.history[1].url, "www.db2URLgoeshere.com")
        self.assertEquals(a.history[1].description, {"db2_id": 12424})

        # invalid nodes should not work
        self.assertRaises(StandardError, StructureNL, self.s, self.hulk,
                          history=[self.invalid_node])

        # too many nodes should not work
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          history=[self.valid_node] * 1000)

    def test_data(self):
        # Structure data is OK due to PMGEncoder/Decoder
        a = StructureNL(self.s, self.hulk, data={"_structure": self.s2})
        self.assertEqual(a.data["_structure"], self.s2,
                         'Data storage is broken')
        self.assertRaises(ValueError, StructureNL, self.s, self.hulk,
                          data={"bad_key": 1})

    def test_eq(self):
        # test basic equals()
        created_at = datetime.datetime.now()
        a = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        b = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        self.assertEqual(a, b, "__eq__() method is broken! false negative")

        # change the created at date, now they are no longer equal
        c = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2])
        self.assertNotEqual(a, c, "__eq__() method is broken! false positive")

        # or try a different structure, those should not be equal
        d = StructureNL(self.s2, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": self.s2},
                        [self.valid_node, self.valid_node2], created_at)
        self.assertNotEqual(a, d, "__eq__() method is broken! false positive")

    def test_to_from_dict(self):
        # no complicated objects in the 'data' or 'nodes' field
        a = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": "string"},
                        [self.valid_node, self.valid_node2])
        b = StructureNL.from_dict(a.to_dict)
        self.assertEqual(a, b)
        # complicated objects in the 'data' and 'nodes' field
        complicated_node = {"name": "complicated node",
                            "url": "www.complicatednodegoeshere.com",
                            "description": {"structure": self.s2}}
        a = StructureNL(self.s, self.hulk, ['test_project'], self.pmg,
                        ['remark1'], {"_my_data": {"structure": self.s2}},
                        [complicated_node, self.valid_node])
        b = StructureNL.from_dict(a.to_dict)
        self.assertEqual(a, b,
                         'to/from dict is broken when object embedding is '
                         'used! Apparently PMGJSONEncoding is broken...')


if __name__ == '__main__':
    unittest.main()
