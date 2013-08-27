"""
This module is used to estimate the cost of various compounds. Costs are taken
from the accompanying CSV file in the format of (formula, cost/kg, name,
reference). For compounds with no cost listed, a Phase Diagram style convex hull
optimization is performed to determine a set of compositions that can be mixed
to give the desired compound with lowest total cost.
"""

from __future__ import division
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
    """
    Extends PDEntry to include a BibTeX reference and include language about
    cost
    """

    def __init__(self, composition, cost, name, reference):
        """
        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition
            cost:
                Cost (per mol, NOT per kg) of the full Composition
            name:
                Optional parameter to name the entry. Defaults to the reduced
                chemical formula.
            reference:
                Reference data as BiBTeX string
        """
        super(CostEntry, self).__init__(composition, cost, name)
        if reference and not is_valid_bibtex(reference):
            raise ValueError(
                "Invalid format for cost reference! Should be BibTeX string.")
        self.reference = reference

    def __repr__(self):
        return "CostEntry : {} with cost = {:.4f}".format(self.composition,
                                                          self.energy)


@singleton
class CostDB(object):
    """
    Singleton object that stores the cost data and provides functions for
    analyzing cost
    """

    def __init__(self):
        """ Implementation of the singleton interface """
        self._chemsys_entries = defaultdict(list)
        filename = os.path.join(os.path.dirname(__file__), "cost_db.csv")
        reader = csv.reader(open(filename, "rb"), quotechar="|")
        for row in reader:
            comp = Composition.from_formula(row[0])
            cost_per_mol = float(row[1]) * comp.weight * \
                AMU_TO_KG * AVOGADROS_CONST
            pde = CostEntry(comp.formula, cost_per_mol, row[2], row[3])
            chemsys = "-".join([el.symbol for el in pde.composition.elements])
            self._chemsys_entries[chemsys].append(pde)

    def get_lowest_decomposition(self, composition, override_known=True):
        """
        Get the decomposition leading to lowest cost

        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition
            override_known:
                If the cost for the Composition is listed explicitly, don't
                minimize cost
        Returns:
            Decomposition as a dict of {Entry: amount}
        """

        # first check if cost is explicitly tabulated for this composition
        if override_known:
            chemsys = "-".join([el.symbol for el in composition.elements])
            for entry in self._chemsys_entries[chemsys]:
                if entry.composition.reduced_formula == \
                        composition.reduced_formula:
                    return {entry: 1}

        # otherwise minimize cost via a convex hull analysis
        entries_list = []
        elements = [e.symbol for e in composition.elements]
        for i in range(len(elements)):
            for combi in itertools.combinations(elements, i + 1):
                chemsys = "-".join(sorted(combi))
                entries_list.extend(self._chemsys_entries[chemsys])

        pd = PhaseDiagram(entries_list)
        return PDAnalyzer(pd).get_decomposition(composition)

    def get_cost_per_mol(self, composition, override_known=True):
        """
        Get best estimate of minimum cost/mol based on known data

        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition
            override_known:
                If the cost for the Composition is listed explicitly, don't
                minimize cost
        Returns:
            float of cost/mol
        """

        decomp = self.get_lowest_decomposition(composition, override_known)
        return sum(k.energy_per_atom * v * composition.num_atoms for k, v in
                   decomp.iteritems())

    def get_cost_per_kg(self, comp, override_known=True):
        """
        Get best estimate of minimum cost/kg based on known data

        Args:
            composition:
                Composition as a pymatgen.core.structure.Composition
            override_known:
                If the cost for the Composition is listed explicitly, don't
                minimize cost
        Returns:
            float of cost/kg
        """
        return self.get_cost_per_mol(comp, override_known) / \
            (comp.weight * AMU_TO_KG * AVOGADROS_CONST)