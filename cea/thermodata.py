"""Extended interface to  the NASA Glenn thermodynamic database

This module is built on top of the `thermoinp` module (which serves to
provide low-level access to the database).

"""
import collections

Interval = collections.namedtuple('Interval', 
								  ['bounds',
								  'coeffs',
								  'integration_consts'])

Species = collections.namedtuple('Species', 
								 ['name',
								  'molar_mass',
								  'heat_formation', 
								  'intervals'])
