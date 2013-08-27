import unittest
from pymatgen import Composition
from pymatgen.analysis.cost import CostDB

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Aug 27, 2013'


class CostTest(unittest.TestCase):
    def setUp(self):
        self.cdb = CostDB()
        self.al = Composition.from_formula('Al')
        self.al2 = Composition.from_formula('Al2')
        self.o = Composition.from_formula('O')
        self.al2o3 = Composition.from_formula('Al2O3')

    def test_molsanity(self):
        # cost per mol doubles as formula unit doubles
        alm = self.cdb.get_cost_per_mol(self.al)
        al2m = self.cdb.get_cost_per_mol(self.al2)
        self.assertAlmostEqual(alm * 2, al2m, 8)

    def test_kgsanity(self):
        # cost per kg constant as formula unit doubles
        al = self.cdb.get_cost_per_kg(self.al)
        al2 = self.cdb.get_cost_per_kg(self.al2)
        self.assertAlmostEqual(al, al2, 8)

    def test_hull(self):
        # cost per mol is additive when using the hull
        # this test assumes no Al-O binary data except for Al2O3
        al = self.cdb.get_cost_per_mol(self.al)
        o = self.cdb.get_cost_per_mol(self.o)
        al2o3 = self.cdb.get_cost_per_mol(self.al2o3, override_known=False)

        self.assertAlmostEqual(al2o3, 2 * al + 3 * o, 8)

    def test_override(self):
        # compound cost is lower than element costs for Al2O3
        # this test assumes no Al-O binary data except for Al2O3
        al = self.cdb.get_cost_per_mol(self.al)
        o = self.cdb.get_cost_per_mol(self.o)
        al2o3 = self.cdb.get_cost_per_mol(self.al2o3)
        self.assertLess(al2o3, 2 * al + 3 * o)



