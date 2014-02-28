import os
import re
from collections import namedtuple

_THERMOINP = os.path.join('cea', 'data', 'thermo.inp')

_Species = namedtuple('ChemicalSpecies',
					  'name, comments, no_intervals, refcode,'
					  'formula, phase, molwt, heat_formation,'
					  'intervals'
					  )

class ChemSpecies(_Species):
	"""Chemical species with encapsulated thermodynamic data."""

	@classmethod
	def from_records(cls, records):
		"""Construct instance from relevant records in thermo.inp"""
		# Each set of records contains 3 subsets
		independent_data, interval_data = records[:2], records[2:]

		# The first record is an identifier
		header = independent_data[0]
		name, comments = header[:18].rstrip(), header[18:].rstrip()

		# Temperature independent attributes live in second record
		subheader = independent_data[1]
		# The number of intervals is an integer in set [1-3]
		no_intervals = int(subheader[1])
		# The refcode is an 8-char string
		refcode = subheader[2:10].strip()
		# The formula is a string of (element, integer) pairs
		formula = subheader[10:51]
		# The phase (0 for gas, non-zero for condensed) is an integer
		phase = int(subheader[51:52])
		# Molecular weight is a decimal
		molwt = float(subheader[52:65].lstrip())
		# Heat of formation is a decimal
		heat_formation = float(subheader[65:80].lstrip())

		# Split the per-interval records up and create a TempInterval
		# object with each one, appending to the intervals list.
		intervals = []
		for records in (interval_data[i:i+3]
						for i in xrange(0, no_intervals * 3, 3)):
			intervals.append(TempInterval.from_records(records))

		return cls(name, comments, no_intervals, refcode, formula,
				   phase, molwt, heat_formation, intervals)

_TempInterval = namedtuple('TempInterval', # type
						   'bounds, exponents, offset, ncoeffs,'
						   'coefficients, constants, a8')

class TempInterval(_TempInterval):
	"""Temperature interval with $C_p/R$ polynomial description.
	
	This is a namedtuple (with modification) representation of a
	temperature interval for a chemical species in the `thermo.inp`
	thermodynamic database. The temperature interval is used for
	curve-fitting purposes allowing evaluation of thermally perfect
	chemical properties at specified temperatures.
	
	"""

	_regex = re.compile(r'\s+')

	@staticmethod
	def _str_to_floats(string):
		float_strings = [string[i:i+16].replace('D','e') # Pythonify
				  		 for i in xrange(0, len(string), 16)]
		return map(float, float_strings)

	@classmethod
	def from_records(cls, records):
		"""Create a temperature interval from `thermo.inp` records

		Database Specification
		----------------------
		
		Interval data is split across three records:
		
		1.	Temperature range
			Number of coefficients for $C_p / R$
			$T$ exponents in empirical equation for $C_p / R$
			${H(298.15) - H(0)}$ J/mol
		2.	Coefficients $a_1 ... a_5$
		3.	Coefficients $a_6, a_7, a_8$
			Integration constants $b_1$ and $b_2$

		Note that here $H$ and $C_p$ refer to molar standard-state
		enthalpy and pressure respectively, as defined in NASA-RP-1311
		Part I, p45. The 8th coefficient, $a_8$, is not treated as a
		coefficient as it's usually blank, or 0.0 if not.
		
		"""

		# Break first record down into constituent parts
		header = records[0] # for clarity
		# Bounds consist of two whitespace-separated decimals
		bounds = tuple(map(float,
						   cls._regex.split(header[1:22].strip())))
		# The number of coefficients is an integer
		ncoeffs = int(header[22:23])
		# The exponents are of the form '-2.0', '1.0' etc, so we need
		# to convert this to a float before integer (which is really
		# what they are, fractional exponents are negative integers).
		exponents = tuple(map(lambda s: int(float(s)),
					    	  cls._regex.split(header[23:63].strip())))
		# The offset {H(298.15) - H(0)} is a decimal
		offset = float(header[65:80])

		# The second and third record are composed of 16-char doubles
		# (with Fortran D-notation).
		coefficients = cls._str_to_floats(records[1])
		coefficients.extend( cls._str_to_floats(records[2][:32]) )
		
		# The 8th coefficient is not used, so treat separately for now
		# TODO: confirm this is never significant, tentatively YAGNI
		try:
			a8 = cls._str_to_floats(records[2][32:48])
		except ValueError:
			a8 = None

		# The remaining doubles are the integration constants b1, b2
		constants = cls._str_to_floats(records[2][48:])

		# Instantiate
		return cls(bounds, exponents, offset, ncoeffs, coefficients,
				   constants, a8)
	
	@property
	def Tmin(self):
		"""The interval's upper temperature limit"""
		return self.bounds[0]

	@property
	def Tmax(self):
		"""The interval's lower temperature limit"""
		return self.bounds[1]
