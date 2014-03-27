"""Interface for NASA Glenn thermodynamic database source file.

This module provides low-level access to the original database format.
The database is broadly split into three categories of chemical
species; gaseous equilibrium products, condensed equilibrium products
and reactants (e.g. 'Air' or other mixtures with characteristic
properties).

Several API functions are included to extract data from the source
file at several levels of decomposition.

  - `read_categories` and `read_species` provide category-keyed
	dictionaries of raw strings and lists of per-species strings
	respectively.
  - `parse` is similar to `read_species` but instead parses the
	per-species strings into a container with field attributes
	(namedtuple).

"""
import re
import os
import collections

def read_categories():
	"""Split the database into categories.
	
	Returns a category-keyed dictionary of string values.

		>>> d = read_categories()
		>>> sorted(d.keys())
		['condensed_products', 'gas_products', 'reactants']
		>>> type(d['reactants'])
		str
	
	"""
	path = os.path.join(os.path.dirname(__file__),
					    'data',
						'thermo.inp')
	keys = 'gas_products', 'condensed_products', 'reactants'
	with open(path, 'r') as f: contents = f.read()
	# Gaseous reactants/products begin with species 'e-'
	# Condensed reactants/products begin with species 'Ag(cr)'
	# Reactants begin with 'Air'
	# We can catch all three categories in one go with a regex:
	pattern = re.compile(r'\n(e-.*)'
		 			      '\n(Ag\(cr\).*)'
			   	    	  '\nEND PRODUCTS.*'
			   		      '\n(Air.*)'
			   		      '\nEND REACTANTS', re.DOTALL)
	return dict(zip(keys, pattern.search(contents).groups()))

def read_species():
	"""Split the database into categorised lists of species.
	
	Returns a category-keyed dictionary of species-separated content. 
	Each category is a list of species dataset strings. The strings
	contain newline-delimited records containing species-specific
	data.

		>>> d = read_species()
		>>> sorted(d.keys())
		['condensed_products', 'gas_products', 'reactants']
		>>> type(d['reactants'])
		list
		>>> type(d['reactants'][0])
		str

	"""
	categories = read_categories()
	# For each category, separate into per-species strings
	pattern = re.compile(r'\n(?=[eA-Z(])')
	return {k:pattern.split(v) for k, v in categories.items()}

def parse():
	"""Parse the database into categorised lists of Species instances.

	Returns a category-keyed dictionary of lists containing Species
	instances.

		>>> d = parse()
		>>> sorted(d.keys())
		['condensed_products', 'gas_products', 'reactants']
		>>> type(d['reactants'])
		list
		>>> sample_species = d['reactants'][0]
		>>> sample_species.name
		'Air'

	"""
	species_categories = read_species()
	return {k:[_parse_species(string.split('\n')) for string in lst]
			for k, lst in species_categories.items()
			}


# --------------------------------------------------------------------
#
# Data structure for the 2002-spec species datasets
#
# --------------------------------------------------------------------

_Species = collections.namedtuple('ChemicalSpecies',
							      ['name',
							       'comments',
							       'nintervals',
							       'refcode',
							       'formula',
							       'phase',
							       'molwt',
								   'h_formation',
								   'h_assigned',
								   'T_reference',
								   'intervals'])

class Species(_Species):
	"""Chemical species metadata and thermodynamic properties.

    	`name` 		  : Species name/ID (usually formula)
    	`comments`    : References and comments
    	`nintervals`  : Number of temperature intervals (0-3)
    	`refcode`     : Reference-date code
    	`formula`     : String array of elements and no. atoms
    	`phase`       : Condensed phases non-zero
    	`molwt`       : molecular weight/molar mass, kg/kmol
    	`h_formation` : heat/enthalpy of formation (nintervals > 0)
    	`h_assigned`  : assigned enthalpy (nintervals == 0)
		`T_reference` : reference temperature (for assigned enthalpy)
		`intervals`   : temperature intervals (nintervals > 0)

    """
	pass


def _parse_species(records):
	# Parse records containing species data.
	# Returns a Species instance.

	# split the records up
	head, body, tail = records[0], records[1], records[2:]

	# Parse the name & comments from the header 
	name = head[:18].rstrip()
	comments = head[18:].rstrip()

	# Parse the non-polynomial data
	nintervals = int(body[1])
	refcode = body[2:10].strip()
	# make the formula a bit more parse-friendly but leave as a string
	# e.g. 'C   1.00O  2.00   0.00   0.00   0.00' -> 'C:1.00 O:2.00'
	formula = ' '.join([
		'{!s}:{!s}'.format(body[i:i+2].strip(), body[i+2:i+8].strip())
		for i in xrange(10, 50, 8)
		]).replace(' :0.00', '')
	phase = int(body[51])
	molwt = float(body[52:65])

	# At this stage, the fields have the potential to vary depending
	# on whether temperature intervals are present or not.
	refenthalpy = float(body[65:])
	if nintervals > 0:
		h_assigned = T_reference = None
		h_formation = refenthalpy
		# each interval is described by three records
		intervals = [_parse_interval(tail[i:i+3])
					 for i in xrange(0, len(tail), 3)]
	else:
		h_formation = intervals = None
		h_assigned = refenthalpy
		T_reference = tail[0].split()[0] # grab the first word

	return Species(name,
				   comments,
				   nintervals,
				   refcode,
				   formula,
				   phase,
				   molwt,
				   h_formation,
				   h_assigned,
				   T_reference,
				   intervals)


# --------------------------------------------------------------------
#
# Data structure for the 2002-spec polynomial form (variable form
# <=8-term poly with 2 integration constants)
#
# --------------------------------------------------------------------

_Interval = collections.namedtuple('TemperatureInterval',
								   ['bounds',
								    'ncoeff',
								    'exponents',
								    'deltah',
								    'coeff',
								    'const'])


class Interval(_Interval):
	"""Specification of polynomial function of temperature.

	This class of object stores data describing a variable form
	polynomial function applicable to a defined temperature interval.
	Fields correspond to the following data:

	  `bounds`    : Interval bounds (Tmin, Tmax)
	  `ncoeff`    : Number of coefficients/terms
	  `exponents` : Exponent magnitudes (len() == ncoeff)
	  `deltah`    : Reference enthalpy value
	  `coeff`     : Coefficients (len() == ncoeff)
	  `const`     : Integration constants
	
	"""
	pass

def _double_array_to_float(string):
	# Parse a string a containing 16-char Fortran-style doubles into
	# a list of floats
	float_strings = [string[i:i+16].replace('D','e') # Pythonify
			  		 for i in xrange(0, len(string), 16)]
	return map(float, float_strings)

def _parse_interval(records):
	# Parse records containing a temperature interval/polynomial spec.
	# This expects records as a list of strings and returns an
	# Interval instance.
	metadata, array1, array2 = records

	# parse metadata string first
	bounds = tuple(float(n) for n in metadata[:22].split())
	ncoeffs = int(metadata[22])
	exponents = tuple(float(n) for n in metadata[23:63].split())
	deltah = float(metadata[65:])

	# parse records containing numerical strings 
	coeffs = _double_array_to_float(array1)
	coeffs.extend(_double_array_to_float(array2[:32]))
	coeffs = tuple(coeffs)
	consts = tuple(_double_array_to_float(array2[48:]))

	return Interval(bounds, 
					ncoeffs,
					exponents,
					deltah,
					coeffs,
					consts)
