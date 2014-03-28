Design Implementation
=====================

These notes provides implementation details for this project.


### Database Interface

The database interface should provide the following features:

  - A source database format parser.
  - Alternative database format(s) (including writer/parser).
  - Basic searching, browsing and creation of data subsets in both the
    original and alternative formats.

Really this deserves two levels of abstraction in my view.

1.	The first is a simple wrapper around the original database. The
	original format is generally respected including redundant data,
	and functionality is limited to simple parsing and basic views on
	the data. This really only needs a single Python module.
2.	The second builds on this, including interfacing with alternative
	data storage/interchange formats, providing functions for
	calculating properties, etc. This will be at least one Python
	module.

The following is a simple representation of the relationships between
system elements:

	alt.db(s) <--> thermodata <--> thermoinp <--> db


#### Source parser

This is a feature of the module `thermoinp`.

The original format is a relatively monolothic formatted text file.
The parser feature should, simply, provide processed contents through
the module's API. This can be done via lightweight data structures
representing the chemical species. One additional element of metadata
worth preserving is the broad categorisation of the database into
gaseous reactants/products, condensed reactants/products and
reactants. 


#### Alternative database formats

This is a feature of the module `thermodata` and builds on the
functionality of the `thermoinp` parser.

The source format has a certain amount of both built-in redundancy
and, due to the benefits of (more) sophisticated data structures,
redundant metadata. In short this means the data is amenable to
representation in a simplified, semantically richer format both
internally within the code and externally in an alternative format.

An attractive option from the perspective of data interchange is XML.
A suitable XML schema is capable of representing semantically rich 
data structures and downsides (such as relatively awkward parsing and
querying) are mitigated by the static nature of the data, relatively
small size and popular XML libraries. 

#### Querying

Querying features can be implemented in both the low-level
(`thermoinp`) and higher-level (`thermodata`) interface modules. The
former allows for generating database subsets based on simple
filtering criteria (name-matching, categories) and the latter offers
opportunity for slightly more advanced operations similar to those
offered by [ThermoBuild][], including generation of tabular data and
plotting capabilities.
