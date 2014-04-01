"""Extended interface to  the NASA Glenn thermodynamic database

This module is built on top of the `thermoinp` module (which serves to
provide low-level access to the database).

"""
from math import log
import collections

import constants as CONST
import thermoinp


Interval = collections.namedtuple('Interval', 
								  ['bounds',
								  'coeffs',
								  'integration_consts'])


class Species(object):
	"""Chemical species"""
	def __init__(self, name):
		try:
			source = db[name]
		except KeyError as e:
			raise e("Unrecognised species: {}".format(name))
		self.name = source
		self.Mr = source.molwt
		self.M = CONST.M * self.Mr
		self.R = CONST.R_CEA / self.M
		self.Hf = source.h_formation
		self.hf = source.h_formation / self.M
		try:
			intervals = tuple(Interval(interval.bounds,
									   interval.coeff, 
									   interval.const)
    	                      for interval in source.intervals
							  )
		except TypeError:
			intervals = None
		self.thermo = Thermo(self, intervals)


class Thermo(object):
	"""Thermodynamic state functions (standard-state, P=100 kPa)."""
	def __init__(self, species, intervals, T=298.15):
		self.species = species
		self.intervals = intervals
		self.bounds = intervals[0].bounds[0], intervals[-1].bounds[1]
		self.T = T
	
	@property
	def T(self):
		"""Temperature, K"""
		return self._T
	@T.setter
	def T(self, T):
		self._T = T

		# Localise variables for repeated access
		Ru = CONST.R_CEA
		R = self.species.R
		a = self.intervals[0].coeffs
		b1, b2 = self.intervals[0].integration_consts

		# Calculate dimensionless values
		Cp_nodim = _dimless_heat_capacity(T, a)
		H_nodim = _dimless_enthalpy(T, a, b1)
		S_nodim = _dimless_entropy(T, a, b2)

		# Assign properties
		self._Cp = Cp_nodim * Ru
		self._cp = Cp_nodim * R
		self._H = H_nodim * Ru * T
		self._h = H_nodim * R * T
		self._S = S_nodim * Ru
		self._s = S_nodim * R


	# Heat capacity properties
	# ----------------------------------------------------------------
	@property
	def Cp(self):
		"""Molar heat capacity at constant pressure, J/mol-K."""
		return self._Cp

	@property
	def cp(self):
		"""Specific heat capacity at constant pressure, J/kg-K."""
		return self._cp

	# Enthalpy properties
	# ----------------------------------------------------------------
	@property
	def H(self):
		"""Molar enthalpy, J/mol"""
		return self._H

	@property
	def h(self):
		"""Specific enthalpy, J/kg"""
		return self._h

	# Entropy properties
	# ----------------------------------------------------------------
	@property
	def S(self):
		"""Molar entropy, J/mol-K"""
		return self._S

	@property
	def s(self):
		"""Specific entropy, J/kg-K"""
		return self._s


def _dimless_heat_capacity(T, a):
	# Returns the dimensionless heat capacity, Cp/R
	# T : Temperature, K
	# a : coefficients, len(a) == 7
	return (  a[0] / T**2 
	        + a[1] / T
	        + a[2]
	        + a[3] * T
	        + a[4] * T**2
	        + a[5] * T**3
	        + a[6] * T**4
	        ) 

def _dimless_enthalpy(T, a, b):
	# Returns the dimensionless enthalpy, H/RT
	# T : Temperature, K
	# a : coefficients, len(a) == 7
	# b : integration constant
	return (- a[0] / T**2
	        + (a[1] * log(T) + b) / T
	        + a[2]
	        + a[3] * T / 2.0 
	        + a[4] * T**2 / 3.0 
	        + a[5] * T**3 / 4.0 
	        + a[6] * T**4 / 5.0
	        )

def _dimless_entropy(T, a, b):
	# Returns the dimensionless entropy, S/R.
	# T : Temperature, K
	# a : coefficients, len(a) == 7
	# b : integration constant
	return (- a[0] / T**2 / 2.0 
	        - a[1] / T 
	        + (a[2] * log(T) + b)
	        + a[3] * T
	        + a[4] * T**2 / 2.0 
	        + a[5] * T**3 / 3.0 
	        + a[6] * T**4 / 4.0
	        )

def _specific_gas_constant(M):
	# Returns the specific gas constant as a function of molar mass
	# M : Molar mass, kg/mol
	return CONST.R_CEA / M


def _load_database():
	# Return the NASA database as a flat dictionary.
	categorised_dict = thermoinp.parse()
	return {species.name:species
		    for species_list in categorised_dict.values()
		    for species in species_list
		    }

db = _load_database()
