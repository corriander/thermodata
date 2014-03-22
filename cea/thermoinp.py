import re
import os

def read():
	"""Read the NASA thermodynamic database into tuple of strings.
	
	Returns (gas_products, condensed_products, reactants)
	
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

def parse():
	"""Return the database into a category-keyed dictionary.
	
	Categories are 'gas', 'condensed' and 'reactant' corresponding to
	gaseous reactant/products, condensed reactant/products and
	reactants (of any state) respectively.
	
	"""
	db = dict.fromkeys(('gas', 'condensed', 'reactant'))
	db['gas'], db['condensed'], db['reactant'] = read()
	# For each category, separate into species then split into records
	pattern = re.compile(r'\n(?=[eA-Z(])')
	for category in db: 
		# Convert the string in a per-species list of strings
		db[category] = pattern.split(db[category])
		for i, species in enumerate(db[category]):
			# convert the species string into a list of records
			species = species.split('\n') 
			db[category][i] = species # TODO parse into a dict
	return db





		
