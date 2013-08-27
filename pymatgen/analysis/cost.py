from collections import defaultdict
import csv
import os
import itertools
from pymatgen import Composition, singleton, AMU_TO_KG
from pymatgen.core.physical_constants import AVOGADROS_CONST
from pymatgen.matproj.snl import is_valid_bibtex
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.pdmaker import PhaseDiagram

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Aug 27, 2013'


class CostEntry(PDEntry):

    def __init__(self, composition, energy, name, reference):
        super(CostEntry, self).__init__(composition, energy, name)
        if reference and not is_valid_bibtex(reference):
            raise ValueError("Invalid format for cost reference! Should be BibTeX string.")
        self.reference = reference

    def __repr__(self):
        return "CostEntry : {} with cost = {:.4f}".format(self.composition,
                                                          self.energy)

@singleton
class CostDB(object):
    """
    """

    def __init__(self):
        """ Implementation of the singleton interface """
        self._chemsys_entries = defaultdict(list)
        self._comp_energies = {}
        filename = os.path.join(os.path.dirname(__file__), "cost_db.csv")
        reader = csv.reader(open(filename, "rb"))
        for row in reader:
            comp = Composition.from_formula(row[0])
            cost_per_mol = float(row[1]) * comp.weight * AMU_TO_KG * AVOGADROS_CONST
            pde = CostEntry(comp.formula, cost_per_mol, row[2], row[3])
            chemsys = "-".join([el.symbol for el in pde.composition.elements])
            self._chemsys_entries[chemsys].append(pde)
            self._comp_energies[comp.reduced_formula] = min(pde.energy_per_atom, self._comp_energies.get(comp.reduced_formula, float('inf')))

            print self._comp_energies

    def get_cost_per_mol(self, composition):
        # first check if cost is explicitly defined for this composition
        if composition.reduced_formula in self._comp_energies:
            return self._comp_energies[composition.reduced_formula] * composition.num_atoms

        # otherwise create a phase diagram to get lowest cost
        entries_list = []
        elements = [e.symbol for e in composition.elements]
        for i in range(len(elements)):
            for combi in itertools.combinations(elements, i + 1):
                chemsys = "-".join(sorted(combi))
                entries_list.extend(self._chemsys_entries[chemsys])

        pd = PhaseDiagram(entries_list)
        pda = PDAnalyzer(pd)

        return sum(k.energy*v*composition.num_atoms for k, v in pda.get_decomposition(composition).iteritems())

    def get_cost_per_kg(self, comp):
        return self.get_cost_per_mol(comp)/(comp.weight * AMU_TO_KG * AVOGADROS_CONST)



if __name__ == "__main__":
    cdb = CostDB()
    print cdb.get_cost_per_mol(Composition.from_formula("LiFeO2"))
    print cdb.get_cost_per_kg(Composition.from_formula("LiFeO2"))