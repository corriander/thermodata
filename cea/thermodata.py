"""Extended interface to  the NASA Glenn thermodynamic database

This module is built on top of the `thermoinp` module (which serves to
provide low-level access to the database).

"""
import sys
from math import log
import collections
from xml.etree import ElementTree as etree

import constants as CONST
import thermoinp


Interval = collections.namedtuple('Interval', 
								  ['bounds',
								  'coeffs',
								  'integration_consts'])


class ChemDB(dict):
	"""Chemical database.

	This class loads the source database in the background on
	isntantiation and then provides a mechanism for selecting species
	from the source to add to the current "view". 
	
		>>> db = ChemDB()
		>>> db
		{}
		>>> db.select(('Air',)) # NOTE: names passed as sequence
		>>> air = db['Air']

	The method will raise an exception if the chemical species does
	not exist in the source:

		>>> db.select(('Adamantium',))
		Traceback (most recent call last):
			...
		Exception: Adamantium is not in the source database.

	The current database entries can be written to `stdout` via the
	`write` method or to a (new) file by specifying the path as an
	argument.

	"""
	def __init__(self):
		self._thermoinp_load()
	
	def select(self, names):
		"""Generate database by a list of species names."""
		for name in names:
			try:
				self.__setitem__(name, self._source_dict[name])
			except KeyError as e:
				msg = "{} is not in the source database.".format(name)
				raise Exception(msg)

	def _thermoinp_load(self):
		# Database loader. Loads the contents of `thermo.inp` into a
		# flat dictionary.
		categ_dict = thermoinp.parse()
		self._source_dict = {species.name:self._map_species(species)
			    			 for species_list in categ_dict.values()
			    			 for species in species_list
			    			 }
	
	def toxml(self):
		"""Represent database contents in XML form."""
		root = etree.Element('chemdb')
		for species_obj in self.values():
			species_obj.toxml(root)

		return root
	
	def write(self, path=None):
		"""Write database in XML format.
		
		  - If file path is unspecified, XML is written to STDOUT.
		  - If file exists, an exception is raised.
		
		"""
		if path is None:
			f = sys.stdout
		else:
			if os.path.isfile(path):
				raise IOError("{} exists.".format(path))
			f = open(path, 'w')

		root = self.toxml()
		_indentxml(root)
		etree.ElementTree(root).write(f,
									  xml_declaration=True,
									  encoding='utf-8',
									  method='xml')
		if f is not sys.stdout:
			f.close()

	def _map_species(self, source):
		# map thermoinp.Species instance data to Species instances.
		try:
			intervals = [self._map_interval(i)
						 for i in source.intervals]
		except TypeError:
			intervals = None
	
		return Species(source.name, source.molwt, source.h_formation,
				   	   intervals)
	
	@staticmethod
	def _map_interval(source):
		# map thermoinp.Interval instance data to Interval instances
		return Interval(source.bounds, source.coeff, source.const)


class Species(object):
	"""Chemical species.
	
	A chemical species is a substance of known composition in a
	defined phase and, if an ion, with a specified charge.
	
	Species can be instantiated directly, but is generally
	instantiated in the database loading during the instantiation of
	ChemDB.
	
	"""
	def __init__(self, name, rel_molar_mass, formation_enthalpy,
				 intervals=None):
		self.name = name
		self.Mr = rel_molar_mass
		self.Hf = formation_enthalpy

		# Derived attributes:
		self.M = CONST.M * self.Mr
		self.R = CONST.R_CEA / self.M
		try:
			self.hf = self.Hf / self.M
		except TypeError:
			# formation enthalpy undefined in source Species
			self.hf = None

		# Attach thermodynamic property model
		if intervals:
			self.thermo = Thermo(self, intervals)
		else:
			self.thermo = None
	
	def toxml(self, parent):
		"""Create an XML representation of the thermodynamic model"""
		attributes = {'name' : self.name}
		node = etree.SubElement(parent, 'species', attributes)
		M = etree.SubElement(node,
							 'molar_mass', 
							 {'units' : 'kg/mol'}
							 )
		M.text = str(self.M)
		R = etree.SubElement(node,
							 'gas_constant',
							 {'units' : 'J/kg-K'}
							 )
		R.text = str(self.R)
		Hf = etree.SubElement(node,
							 'formation_enthalpy',
							 {'units' : 'J/mol'}
							 )
		Hf.text = str(self.Hf)
		self.thermo.toxml(node)


class Thermo(object):
	"""Thermodynamic state functions (standard-state, P=100 kPa).
	
	Temperature, T, is used as the free variable here. On setting the
	temperature property state functions are evaluated and stored in
	attributes. In the event no intervals are provided, this
	evaluation process does not happen (the necessary data is not
	available).

	The following properties are available for standard-state
	conditions (specified temperature and standard pressure, 
	P =	100 kPa):

	  - T      : Temperature
	  - Cp, cp : Standard-state heat capacity at constant pressure 
	  - H, h   : Standard-state enthalpy
	  - S, s   : Standard-state entropy

	Note that upper-case and lower-case properties are in units of
	amount-of-substance (/mol) and mass (/kg) respectively.
	
	Like Species, Thermo can be instantiated directly but is generally
	handled during the ChemDB database loading.
		
	"""
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
		self._select_interval(T)

		# Localise variables for repeated access
		Ru = CONST.R_CEA
		R = self.species.R
		if self.interval:
			a = self.interval.coeffs
			b1, b2 = self.interval.integration_consts

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

	def _select_interval(self, T):
		# Select the appropriate interval for the current temperature
		if self.bounds[1] < T  < self.bounds[0]:
			raise ValueError("Temperature out of bounds")
		self.interval = None
		for interval in self.intervals:
			if T < interval.bounds[1]:
				self.interval = interval
				break		
	
	def toxml(self, parent):
		"""Create an XML representation of the thermodynamic model"""
		attributes = {'Tmin' : str(self.bounds[0]),
					 'Tmax' : str(self.bounds[1])}
		node = etree.SubElement(parent, 'thermo', attributes)
		for interval in self.intervals:
			attributes = {'Tmin' : str(interval.bounds[0]),
						  'Tmax' : str(interval.bounds[1])
						  }
			subnode = etree.SubElement(node, 'interval', attributes)
			coeffs = etree.SubElement(subnode, 'coefficients')
			coeffs.text = '{!s}'.format(interval.coeffs)
			consts = etree.SubElement(subnode, 'integ_constants')
			consts.text = '{!s}'.format(interval.integration_consts)


class Table(object):
	"""Tabulated data.
	
	Generates tabulated state function values for a specified
	temperature range and chemical species.

		>>> db = ChemDB()
		>>> db.select(('Air',))
		>>> air = db['Air']
		>>> T_range = (200, 500, 2000)
		>>> table = Table(T_range, air)
		>>> str(table)
		'Property table: T = 200-2000 K, 3 intervals (moles)'
		>>> print table.formatted()
		         T        Cp    H-H298         S         H
		         K   J/mol-K    kJ/mol   J/mol-K    kJ/mol
		--------------------------------------------------
		  200         29.034    -2.852   187.221    -2.978
		  500         29.821     5.932   214.001     5.807
		  2000        36.216    56.595   259.764    56.470
	
	"""
	def __init__(self, temperature_range, species):
		self.Trange = temperature_range
		self.species = species
		self._tabulate()
	

	def _tabulate(self):
		# Produced tabulated data.
		header = ('T', 'Cp', 'H-H298', 'S', 'H')
		units = ('K', 'J/mol-K', 'kJ/mol', 'J/mol-K', 'kJ/mol')
		body = []
		species = self.species
		for T in self.Trange:
			species.thermo.T = T
			Cp = species.thermo.Cp
			H = species.thermo.H / 1000
			S = species.thermo.S
			HH298 = (H - species.Hf / 1000)
			row = T, Cp, HH298, S, H
			body.append(row)
		self.header, self.units, self.body = header, units, body
	
	def __str__(self):
		# print a table summary
		Trange = self.Trange
		T_str = 'T = {}-{} K, {} intervals'.format(Trange[0],
												   Trange[-1],
												   len(Trange)
												   )
		units = '(moles)'
		return 'Property table: {} {}'.format(T_str, units)

	def formatted(self):
		"""Format the table for printing/writing to file."""
		spec = '{:>10}' # right-aligned, column width 9
		header = ''.join(spec.format(field) for field in self.header)
		units = ''.join(spec.format(units) for units in self.units)
		fspec = '{:>10.3f}' 
		table = [header, units, '-'*len(header)]
		for row in self.body:
			row = ('  {:<8}'.format(row[0]),
				   ''.join(fspec.format(value) for value in row[1:])
				   )
			table.append(''.join(row))
		return '\n'.join(table)


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

def _indentxml(elem, level=0):
	# Indent XML string representation of elements;
	# http://effbot.org/zone/element-lib.htm#prettyprint
	indent = "    "
	i = "\n" +level*indent
	if len(elem):
		if not elem.text or not elem.text.strip():
			elem.text = i + indent
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
		for elem in elem:
			_indentxml(elem, level+1)
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail = i

if __name__ == '__main__':
	import doctest
	doctest.testmod()
