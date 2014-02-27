import os

_THERMOINP = os.path.join('cea', 'data', 'thermo.inp')

class ChemSpecies(object):

	# Figure out the best way to deal with these attributes
	comments = None # Optional, raw values might be amenable to processing
	refcode = None # Internal reference code, optional esp. without context
	intervals = None # List of TempInterval objects...
	no_intervals = None # Probably derive from intervals unless specified

	"""Chemical species with encapsulated thermodynamic data."""
	def __init__(self, name, formula):
		self.name = name # UID
		self.formula = formula 	# Just the raw string for now
								# It does make the chemical
								# unambiguous though, so probably keep
								# it as required.

	@classmethod
	def from_records(cls, records):
		"""Construct instance from relevant records in thermo.inp"""
		pass

class TempInterval(object):
	"""Temperature interval with $C_p/R$ polynomial description."""
	def __init__(self, 
			T_bounds, offset, coefficients, constants, poly_spec)
		self.min, self.max = T_bounds # Temperature interval bounds
		self.offset = offset # [H_s(298.15) - H_s(0)]
		self.coefficients = coefficients # coefficients, a_1...a_n
		self.constants = constants  # integration constants b1, b2
		self.poly_spec = poly_spec	# raw string for now

	@classmethod
	def from_records(cls, records):
		"""Create a temperature interval from thermo.inp records"""
		pass
