import unittest
import collections
from ..thermodata import Interval, Species, Thermo, ChemDB, Table

class TestSpecies(unittest.TestCase):
	"""Test Species instantiated w/ and w/o formation_enthalpy."""
	def setUp(self):
		self.species = Species('Propane', 44.09562, -104680.0)
		self.species_with_Href = Species('Species', 2.00, None)

	def test_M(self):
		"""Test that M (molar mass) is derived from Mr correctly."""
		self.assertEqual(self.species.M, 0.0441, delta=0.0001)
		self.assertEqual(self.species_with_Href.M, 0.002, delta=0.0001)
	
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

	def test_assign_subbound_T(self):
		"""Tests an exception is raised when T << T_min,species"""
		self.assertRaises(ValueError,
						  lambda t: self.thermo.__setattr__('T', t),
						  0.1)

	def test_assign_superbound_T(self):
		"""Tests a ValueError is raised when T >> T_max,species"""
		self.assertRaises(ValueError,
						  lambda t: self.thermo.__setattr__('T', t),
						  30000.0)
	
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
		db = ChemDB().select(species)
		species = [db[s] for s in species]

		# Generate the tables
		self.table = [
			Table((200, 298.15, 500, 1000, 3000, 6000, 10000, 20000),
				  species[0]),
			Table((200, 298.15, 500, 1000, 3000, 6000), species[1]),
			Table((100, 298.15, 400), species[2]),
			Table((200, 298.15, 500, 1000, 3000, 6000), species[3])
			]

	def fetch_table(self, species):
		"""Utility function to fetch the table being tested."""
		fname = 'table{}.txt'.format(species)
		path = os.path.join(os.path.dirname(__file__), 'data', fname)
		with open(path, 'r') as f: return f.read()

	def test_species_0(self):
		"""Test a gaseous reactant/product with 3 regular intervals"""
		reference_table = self.fetch_table(table[0].species.name)
		self.assertEqual(self.table[0].formatted(), reference_table)

	def test_species_1(self):
		"""Test a gaseous reactant/product with 2 regular intervals"""
		reference_table = self.fetch_table(table[1].species.name)
		self.assertEqual(self.table[1].formatted(), reference_table)

	def test_species_2(self):
		"""Test a condensed species with non-standard breakpoints."""
		reference_table = self.fetch_table(table[2].species.name)
		self.assertEqual(self.table[2].formatted(), reference_table)

	def test_species_3(self):
		"""Test a regular reactant"""
		reference_table = self.fetch_table(table[3].species.name)
		self.assertEqual(self.table[3].formatted(), reference_table)


class TestChemDB(unittest.TestCase):
	# TODO test select for single and then multiple species. Check
	# identity (this catches Species parsing integrity!) 
	pass


if __name__ == '__main__':
	unittest.main()
