"""Extended interface to  the NASA Glenn thermodynamic database

This module is built on top of the `thermoinp` module (which serves to
provide low-level access to the database).

"""
from math import log
import collections

import constants as CONST

Interval = collections.namedtuple('Interval', 
								  ['bounds',
								  'coeffs',
								  'integration_consts'])

Species = collections.namedtuple('Species', 
								 ['name',
								  'molar_mass',
								  'heat_formation', 
								  'intervals'])


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
	return (- coeffs[0] / T**2 / 2.0 
	        - coeffs[1] / T 
	        + (coeffs[2] * log(T) + consts[1])
	        + coeffs[3] * T
	        + coeffs[4] * T**2 / 2.0 
	        + coeffs[5] * T**3 / 3.0 
	        + coeffs[6] * T**4 / 4.0
	        )

def _specific_gas_constant(M):
	# Returns the specific gas constant as a function of molar mass
	# M : Molar mass, kg/mol
	return CONST.R_CEA / M
