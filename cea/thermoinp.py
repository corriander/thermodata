import re
import os

def read_categories():
	"""Split the database into categories.
	
	Returns three strings (gaseous products/reactants, condensed
	products/reactants and mixed-state reactants).

		>>> gases, condensed, reactants = read_categories()
		>>> type(gases)
		<type 'str'>
	
	"""
	path = os.path.join(os.path.dirname(__file__),
					    'data',
						'thermo.inp')
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
	return pattern.search(contents).groups()

def read_species():
	"""Split the database into categories and species.
	
	Returns three lists (gaseous products/reactants, condensed
	products/reactants and mixed-state reactants) containing
	per-species strings.

		>>> gases, condensed, reactants = read_species()
		>>> type(gases)
		<type 'list'>
		>>> type(gases[0])
		<type 'str'>

	"""
	categories = read()
	# For each category, separate into per-species strings
	pattern = re.compile(r'\n(?=[eA-Z(])')
	return map(pattern.split, categories) 
