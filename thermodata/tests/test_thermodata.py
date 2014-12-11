import unittest
import os
import re
import collections
from thermodata.thermodata import Interval, Species, Thermo, ChemDB, Table
from thermodata.thermodata import thermoinp

class TestSpecies(unittest.TestCase):
    """Test Species instantiated w/ and w/o formation_enthalpy."""
    def setUp(self):
        self.species = Species('Propane', 44.09562, -104680.0)
        self.species_with_Href = Species('Species', 2.00, None)

    def test_M(self):
        """Test that M (molar mass) is derived from Mr correctly."""
        self.assertAlmostEqual(self.species.M, 0.0441, delta=0.0001)
        self.assertAlmostEqual(self.species_with_Href.M, 0.002,
                               delta=0.0001)

    def test_R(self):
        """Test that R (spec. gas constant) is derived correctly."""
        self.assertAlmostEqual(self.species.R, 188.6, delta=0.1)
        self.assertAlmostEqual(self.species_with_Href.R, 4157.3,
                               delta=0.1)

    def test_hf(self):
        """Test that the specific enthalpy of formation is derived."""
        self.assertAlmostEqual(self.species.hf, -2.374e6, delta=1e3)
        self.assertIs(self.species_with_Href.hf, None)

    def test_toxml(self):
        """Test xml generated properly."""
        self.assertTrue(False, msg="TBC")


class TestThermo(unittest.TestCase):
    """Test basic functionality of Thermo."""
    def setUp(self):
        # The only information Thermo gets from its associated species
        # is the gas constant currently. Mock up species object.
        MockSpecies = collections.namedtuple('MockSpecies', 'R')
        # Propane as a reference species.
        self.species = MockSpecies(188.556)
        self.intervals = [
            Interval((200.0, 1000.0),
                     (-2.433144337e+05,
                       4.656270810e+03,
                      -2.939466091e+01,
                       1.188952745e-01,
                      -1.376308269e-04,
                       8.814823910e-08,
                      -2.342987994e-11),
                     (-3.540335270e+04,
                       1.841749277e+02)
                     ),
            Interval((1000.0, 6000.0),
                     ( 6.420731680e+06,
                      -2.659791134e+04,
                       4.534356840e+01,
                      -5.020663920e-03,
                       9.471216940e-07,
                      -9.575405230e-11,
                       4.009672880e-15),
                     ( 1.455582459e+05,
                      -2.818374734e+02)
                     )]
        self.thermo = Thermo(self.species, self.intervals)

    def test_bounds (self):
        """Test that the absolute bounds to the data are determined"""
        self.assertEqual(self.thermo.bounds, (200.0, 6000.0))

    def test_default_T(self):
        """Tests that the default temperature is assigned on init"""
        self.assertEqual(self.thermo.T, 298.15)

    def test_assign_negative_T(self):
        """Tests an exception is raised when T < 0.

        T < 0 K is, by definition, invalid.

        """
        self.assertRaises(ValueError,
                          lambda t: self.thermo.__setattr__('T', t),
                          -0.1)

    def test_assign_zero_T(self):
        """Tests an exception is raised when T == 0.

        T = 0 K can potentially cause singularity issues.

        """
        self.assertRaises(ValueError,
                          lambda t: self.thermo.__setattr__('T', t),
                          0.0)
    # TODO: This functionality is currently on hold, default
    # temperature may not be suitable for some species temperature
    # bounds.
    #def test_assign_subbound_T(self):
    #	"""Tests an exception is raised when T << T_min,species"""
    #	self.assertRaises(ValueError,
    #					  lambda t: self.thermo.__setattr__('T', t),
    #					  0.1)

    #def test_assign_superbound_T(self):
    #	"""Tests a ValueError is raised when T >> T_max,species"""
    #	self.assertRaises(ValueError,
    #					  lambda t: self.thermo.__setattr__('T', t),
    #					  30000.0)

    def test_assign_reference_T(self):
        """Tests the reference properties."""
        self.assertEqual(self.thermo.T, 298.15)
        self.assertEqual(self.thermo.interval, self.intervals[0])
        self.assertAlmostEqual(self.thermo.Cp, 7.359e1, delta=1e-2)
        self.assertAlmostEqual(self.thermo.cp, 1.669e3, delta=1)
        self.assertAlmostEqual(self.thermo.H, -1.047e5, delta=1e2)
        self.assertAlmostEqual(self.thermo.h, -2.374e6, delta=1e3)
        self.assertAlmostEqual(self.thermo.S, 2.703e2, delta=1e-1)
        self.assertAlmostEqual(self.thermo.s, 6.130e3, delta=1)
        # TODO: (external) use a (latex?) table to present the values

    def test_interval_a(self):
        """Tests the interval attribute is updated correctly."""
        self.thermo.T = 350.
        self.assertEqual(self.thermo.interval, self.intervals[0])
        self.assertAlmostEqual(self.thermo.Cp, 8.401e1, delta=1e-2)
        self.assertAlmostEqual(self.thermo.cp, 1.905e3, delta=1)
        self.assertAlmostEqual(self.thermo.H, -1.006e5, delta=1e2)
        self.assertAlmostEqual(self.thermo.h, -2.281e6, delta=1e3)
        self.assertAlmostEqual(self.thermo.S, 2.829e2, delta=1e-1)
        self.assertAlmostEqual(self.thermo.s, 6.416e3, delta=1)

    def test_interval_b(self):
        """Tests the interval attribute is updated correctly."""
        self.thermo.T = 1100.
        self.assertEqual(self.thermo.interval, self.intervals[1])
        self.assertAlmostEqual(self.thermo.Cp, 1.827e2, delta=1e-1)
        self.assertAlmostEqual(self.thermo.cp, 4.143e3, delta=1)
        self.assertAlmostEqual(self.thermo.H, 5.664e3, delta=1)
        self.assertAlmostEqual(self.thermo.h, 1.284e5, delta=1e2)
        self.assertAlmostEqual(self.thermo.S, 4.344e2, delta=1e-1)
        self.assertAlmostEqual(self.thermo.s, 9.851e3, delta=1)


class TestTable(unittest.TestCase):
    """Tests a range of species properties in tabular form."""
    # NOTE: Species with a single datum are not handled by ThermoBuild
    def setUp(self):
        # Define the test species
        species = ('CO2', 'C3H8', 'In(cr)', 'Air')
        db = ChemDB()
        db.select(species)
        species = [db[s] for s in species]

        # Generate the tables
        self.table = [
            Table((200, 298.15, 500, 1000, 3000, 6000, 10000, 20000),
                  species[0]),
            Table((200, 298.15, 500, 1000, 3000, 6000), species[1]),
            Table((100, 298.15, 400), species[2]),
            Table((200, 298.15, 500, 1000, 3000, 6000), species[3])
            ]
        self.spaces = re.compile(r'[ \t]+')

    def fetch_table(self, species):
        """Utility function to fetch the table being tested."""
        fname = 'table{}.txt'.format(species)
        path = os.path.join(os.path.dirname(__file__), 'data', fname)
        with open(path, 'r') as f:
            return self.spaces.sub('  ', f.read().rstrip('\n'))

    def normalise_native_table(self, table):
        """Helper function to normalise the format of tables."""
        table = '\n'.join(table.formatted().split('\n')[3:])
        return self.spaces.sub('  ', table)

    @staticmethod
    def table_row_comparison(table, reference_table, msg_prefix=''):
        """Diagnostic method. Print row-by-row comparison."""
        table = table.split('\n')
        reference_table = reference_table.split('\n')

        # Might as well throw an error if table structure
        # fundamentally different. Should be an edge case, but this
        # can be improved to provide further insight if necessary.
        if len(table) != len(reference_table):
            msg_base = "Tables have different num. rows"
            if msg_prefix:
                errmsg = "{}: {}".format(msg_prefix, msg_base)
            else:
                errmsg = msg_base
            raise ValueError(errmsg)

        # return a string containing row mismatches.
        lsofstr = ['',
                   'Erroneous rows ({}):'.format(msg_prefix),
                   '']
        for i, row in enumerate(table):
            ref_row = reference_table[i]
            if row != ref_row:
                lsofstr.extend([row, ref_row, '\n'])
        return '\n'.join(lsofstr)


    def test_species_0(self):
        """Test a gaseous reactant/product with 3 regular intervals"""
        reference_table = self.fetch_table(self.table[0].species.name)
        table = self.normalise_native_table(self.table[0])
        failmsg = self.table_row_comparison(table,
                                            reference_table,
                                            self.table[0].species.name
                                            )
        self.assertEqual(table, reference_table, failmsg)

    def test_species_1(self):
        """Test a gaseous reactant/product with 2 regular intervals"""
        reference_table = self.fetch_table(self.table[1].species.name)
        table = self.normalise_native_table(self.table[1])
        failmsg = self.table_row_comparison(table,
                                            reference_table,
                                            self.table[1].species.name
                                            )
        self.assertEqual(table, reference_table, failmsg)

    def test_species_2(self):
        """Test a condensed species with non-standard breakpoints."""
        reference_table = self.fetch_table(self.table[2].species.name)
        table = self.normalise_native_table(self.table[2])
        failmsg = self.table_row_comparison(table,
                                            reference_table,
                                            self.table[2].species.name
                                            )
        self.assertEqual(table, reference_table, failmsg)

    def test_species_3(self):
        """Test a regular reactant"""
        reference_table = self.fetch_table(self.table[3].species.name)
        table = self.normalise_native_table(self.table[3])
        failmsg = self.table_row_comparison(table,
                                            reference_table,
                                            self.table[3].species.name
                                            )
        self.assertEqual(table, reference_table, failmsg)


class TestChemDB(unittest.TestCase):
    """Test chemical database instantiation."""
    def setUp(self):
        self.db = ChemDB()
        self.methanol = Species('CH3OH(L)', 32.04186, -238910.000,
                intervals=[
                    Interval((175.610, 390.0),
                             (-1.302004763e6, 3.166984180e4,
                              -3.031242152e2, 1.602231130,
                              -4.594507340e-3, 6.990178310e-6,
                              -4.207388950e-9),
                             (-1.656168201e5, 1.514346642e3)
                             )])

    def test_select_all(self):
        """Test all species can be selected."""
        self.db.select()
        self.assertEqual(len(self.db), 2074)

    ## TODO: Feature not implemented yet.
    #def test_select_category(self):
    #	"""Test species can be selected per-category."""
    #	self.db.select('gas_products')
    #	self.assertEqual(len(self.db), 1300)
    #	self.db.select('condensed_products')
    #	self.assertEqual(len(self.db), 2000)
    #	self.db.select('reactants')
    #	self.assertEqual(len(self.db), 2060)

    def test_single_select(self):
        """Test a single species can be selected."""
        self.db.select('CH3OH(L)')

        self.assertEqual(len(self.db), 1)
        self.assertEqual(self.db['CH3OH(L)'], self.methanol)

    def test_multi_select(self):
        """Test multiple species can be selected."""
        species = ('CH3OH(L)', 'KCL', 'Ag(cr)')
        self.db.select(species)

        # Generate extra species for comparison.
        silver_cryst = Species('Ag(cr)', 107.8682, 0.000, intervals=[
            Interval((200.0, 1235.080),
                     (-7.099236470e4, 7.254788020e2, 1.066518380e-1,
                       5.529541550e-3, -4.425590850e-6,
                       2.091668120e-9, -3.888924460e-13),
                     (-4.614014260e3, 5.074216040)
                     )])
        potassium_chloride = Species('KCL', 74.5513, -214574.881,
                intervals=[
                    Interval((200.0, 1000.0),
                             (9.058351510e3, -2.456801212e2,
                              5.680696190, -2.900127425e-3,
                              4.130983060e-6, -2.907340629e-9,
                              8.223850870e-13),
                             (-2.597304884e4, -3.677976854)
                             ),
                    Interval((1000.0, 6000.0),
                             (-2.122945722e5, 9.346158900e2,
                              2.866264958, 1.468386693e-3,
                              -5.834260780e-7, 1.255777709e-10,
                              -9.150148000e-15),
                             (-3.273787640e4, 1.401864636e1)
                             )])

        self.assertEqual(len(self.db), 3)
        self.assertEqual(self.db['CH3OH(L)'], self.methanol)
        self.assertEqual(self.db['KCL'], potassium_chloride)
        self.assertEqual(self.db['Ag(cr)'], silver_cryst)


if __name__ == '__main__':
    unittest.main()
