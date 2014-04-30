import unittest

class TestSpecies(unittest.TestCase):
	"""Test Species instantiated w/ and w/o formation_enthalpy."""
	def setUp(self):
		self.species = Species('Propane', 44.09562, -104680.0)
		self.species_with_Href = Species('Species', 2.00, None)

	def test_M(self):
		"""Test that M (molar mass) is derived from Mr correctly."""
		self.assertEqual(self.species.M, 0.04409562)
		self.assertEqual(self.species_with_Href.M, 0.02)
	
	def test_R(self):
		"""Test that R (spec. gas constant) is derived correctly."""
		self.assertAlmostEqual(self.species.R, 188.56)
		self.assertAlmostEqual(self.species_with_Href.R, 415.73)
	
	def test_hf(self):
		"""Test that the specific enthalpy of formation is derived."""
		self.assertAlmostEqual(self.species.hf, -2.374e6)
		self.assertAlmostEqual(self.species_with_Href.hf, None)

	def test_toxml(self):
		"""Test xml generated properly."""
		self.assertTrue(False, msg="TBC")


if __name__ == '__main__':
	unittest.main()
