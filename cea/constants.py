"""Physical constants.

This module provides the following externally defined physical
constants:

	R     : Molar gas constant (sometimes universal or ideal), J/mol-K
	R_CEA : Molar gas constant defined by Gordon and McBride
	M	  : Molar mass constant, kg/mol 

Note:

The molar gas constant, R, is given by the CODATA 2010 
recommendations. However, CEA employs a slightly different value for
the dimensionless form of the (molar) heat capacity at constant
pressure (the molar enthalpy and entropy functions are integrals of
this):

	Cp/R = sum(xi*T**yj) (i=1:7, j=-2:4)

where R = 8.314510 J/mol-K. In order to consistently evaluate the
molar specific heat capacity (and molar enthalpy and entropy) this
value is required.

References
----------

 - Gordon and McBride, "Computer Program for Calculating Chemical
   Equilibrium with Applications", 1996
 - Mohr, Taylor and Newell, "CODATA Recommended Values of the
   Fundamental Physical Constants: 2010", 2012

"""
import os
from xml.etree import ElementTree as ET

def fetch_value(name):
	"""Return physical constant value from the XML tree"""
	tree = ET.parse(os.path.join(os.path.dirname(__file__), 
								 'data',
								 'constants.xml'
								 )
					)
	root = tree.getroot()
	xpath = "./PhysicalConstant[@name='{}']/value".format(name)
	text = root.find(xpath).text
	return float(text)

M = fetch_value('molar mass constant') 			# kg/mol
R = fetch_value('molar gas constant') 			# J/mol-K
R_CEA = fetch_value('cea molar gas constant')	# J/mol-K
