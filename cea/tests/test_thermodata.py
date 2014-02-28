import unittest
from ..thermodata import TempInterval

sample_species = """
Ar                Ref-Elm. Moore,1971. Gordon,1999.                             
 3 g 3/98 AR  1.00    0.00    0.00    0.00    0.00 0   39.9480000          0.000
    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
 0.000000000D+00 0.000000000D+00 2.500000000D+00 0.000000000D+00 0.000000000D+00
 0.000000000D+00 0.000000000D+00                -7.453750000D+02 4.379674910D+00
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
 2.010538475D+01-5.992661070D-02 2.500069401D+00-3.992141160D-08 1.205272140D-11
-1.819015576D-15 1.078576636D-19                -7.449939610D+02 4.379180110D+00
   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
-9.951265080D+08 6.458887260D+05-1.675894697D+02 2.319933363D-02-1.721080911D-06
 6.531938460D-11-9.740147729D-16                -5.078300340D+06 1.465298484D+03
""".strip()
# Note that the D notation in Fortran denotes double-precision.

expected_species_attributes = {
		'name' : "Ar",
		'comments' : "Ref-Elm. Moore,1971. Gordon,1999.",
		'no_intervals' : 3,
		'refcode' : "g 3/98",
		# TODO: Formula representation is in the works.
		'formula' : "AR  1.00    0.00    0.00    0.00    0.00 ",
		'MW' : 39.9480000,
		'H_f' : 0.000
		}

expected_interval_properties = [
		{
			'min' : 200.000,
			'max' : 1000.000,
			'ncoeffs' : 7,
			'exponents' : (-2, -1, 0, 1, 2, 3, 4, 0),
			'offset' : 6197.428,
			'coefficients' : (
				0.000000000e+00,
				0.000000000e+00,
				2.500000000e+00,
				0.000000000e+00,
				0.000000000e+00,
				0.000000000e+00,
				0.000000000e+00
				),
			'constants' : (
				-7.453750000e+02,
				4.379674910e+00
				),
			'a8' : None,
		},
		{
			'min' : 1000.000,
			'max' : 6000.000,
			'ncoeffs' : 7,
			'exponents' : (-2, -1, 0, 1, 2, 3, 4, 0),
			'offset' : 6197.428,
			'coefficients' : (
				 2.010538475e+01,
				 -5.992661070e-02,
				 2.500069401e+00,
				 -3.992141160e-08,
				 1.205272140e-11,
				 -1.819015576e-15,
				 1.078576636e-19
				 ),
			'constants' : (
				-7.449939610e+02,
				4.379180110e+00
				),
			'a8' : None,
		},
		{
			'min' : 6000.000,
			'max' : 20000.000,
			'ncoeffs' : 7,
			'exponents' : (-2, -1, 0, 1, 2, 3, 4, 0),
			'offset' : 6197.428,
			'coefficients' : (
				-9.951265080e+08,
				6.458887260e+05,
				-1.675894697e+02,
				2.319933363e-02,
				-1.721080911e-06,
				6.531938460e-11,
				-9.740147729e-16,
				),
			'constants' : (
				-5.078300340e+06,
				1.465298484e+03,
				),
			'a8' : None,
		}]

class TestTempInterval(unittest.TestCase):
	def setUp(self):
		self.records = sample_species.split('\n')[2:]
		self.intervals = (
				TempInterval.from_records(self.records[:3]),
				TempInterval.from_records(self.records[3:6]),
				TempInterval.from_records(self.records[6:9])
				)
	
	def test___str_to_floats(self):
		# Use the string starting '  2.010538475D+01-5.99...'
		floats = TempInterval._str_to_floats(self.records[4])
		self.assertItemsEqual(
				floats,
				expected_interval_properties[1]['coefficients'][:5]
				)
	
	def test_bounds(self):
		self.assertEqual(
				self.intervals[1].bounds,
				(1000.000, 6000.000)
				)
	
	def test_min(self):
		self.assertEqual(self.intervals[1].Tmin, 1000.000)
	
	def test_max(self):
		self.assertEqual(self.intervals[1].Tmax, 6000.000)
	
	def test_ncoeffs(self):
		self.assertEqual(self.intervals[1].ncoeffs, 7)
	
	def test_exponents(self):
		self.assertEqual(
				self.intervals[1].exponents,
				expected_interval_properties[1]['exponents']
				)
	
	def test_offset(self):
		self.assertEqual(
				self.intervals[1].offset,
				expected_interval_properties[1]['offset']
				)
	
	def test_coefficients(self):
		self.assertItemsEqual(
				self.intervals[1].coefficients,
				expected_interval_properties[1]['coefficients']
				)
	
	def test_constants(self):
		self.assertItemsEqual(
				self.intervals[1].constants,
				expected_interval_properties[1]['constants']
				)

if __name__ == '__main__':
	unittest.main()
