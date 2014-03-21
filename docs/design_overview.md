Design Overview
===============

This Python code is intended primarily to provide an interface to the
NASA Glenn thermodynamic database distributed with the NASA Glenn
program [Chemical Equilibrium with Applications (CEA)][CEA] by Sanford
Gordon and Bonnie J. McBride. A secondary goal is to provide a Python
interface to CEA itself. This document summarises the rationale and
high-level design approach to achieving this.

Rationale
---------

This code is an extension of a utility I wrote as part of my research
a few years ago in Scilab. Originally I was using CEA to perform
calculate equilibrium chemical composition but required access to a
small subset of the database to perform external calculations. This
started off as a copy-paste type implementation and involved into a
manual method of searching the database for particular chemical
species and retrieving the relevant records. I'd always intended to
extend this to provide more programmatic, less error-prone access but
never had the time to implement it. This is a realisation of that
intent.

The database is an ASCII plain text file with a well-documented format
(see <database.md> for an overview of this, or references at the [CEA
portal][CEA] for full documentation). As such, there is nothing wrong
with that format and it's widely used. However, browsing, searching or
retrieving subsets of the database for specific purposes is limited to
[ThermoBuild][] which is not general-purpose enough for integration
into a custom program. This is an attempt to provide an alternative
data format, interface and functionality to a rich and very useful
dataset.

Design Specification
--------------------

The following is a loose design specification oriented around the two
goals.

### Database Interface

There are several components to this:

  - A source database parser
  - An alternative database format (e.g. XML)
  - Alternative database parser
  - Searching/browsing functionality
  - Methods for creating predefined and custom subsets.
  - Validation, validation, validation.

The database parser should conform to the full database
specifications. The alternative format should probably strip redundant
information from the database, however. XML is a good choice for
inter-application data transfer (especially with XML-parsing libraries
being commonplace in 2014) but performance should be investigated.
There is also scope for a second alternative format (e.g. separating
gaseous species from condensed, products from reactants) though this
ties into predefined subsets. Obviously a parser is required for this
alternative format.

Browsing/searching involves creating suitable species representations
and methods for searching species name and comment data. In this UI
should be methods for sorting/preparing subsets of the data.

Finally, ensuring the integrity of the parsed data is important. Tests
should be built in to provide confidence about this integrity.

### CEA Interface

Simply, the CEA interface needs to provide a mechanism for preparing
input files, reading output files and calling the CEA program
efficiently and cleanly.

There are different ways in which the CEA program can be called:

  - The simplest implementation is a wrapper around the compiled
    executable.
  - The next level is a source modification altering the requirements
    for preparing input files and parsing output files (e.g.  reading
	from `STDIN` or `STDOUT`).
  - Finally, a tight integration with Python might involve direct
    invocation of the FORTRAN subroutines, though this is somewhat new
	ground for me.

Implementation
--------------

This section refers back to the specification details previously
outlined for each component of the code.

### Database Interface

#### Source parser

The role of the source parser is to:

1.	Read the data file
2.	Respect and maintain granular differences in content (chemical
    products, chemical reactants, whether the species is a gas or
	condensed).
3.	Split the database into per-species datasets and convert these
    datasets into python class instances representing each.

The first two roles are suitably covered by two functions, one to read
the significant components of the data file, the other to split the
returned data into reactants/products and reactants only. The database
is broadly split into

  - Gaseous chemical reactants/products
  - Condensed chemical reactants/products
  - Chemical reactants (that are not products, e.g. jet fuel)

Chemical reactants are a comparatively small data subset and may also
be in gaseous or condensed forms. Whilst there is a point to
distinguishing between the gaseous and condensed products, there is
less point splitting this small dataset the same way (the original
database doesn't). 

As for the final role of the parser, splitting the records into
per-species datasets is a little trickier: Species datasets have a
variable number of records which can be determined either through
inspection of the dataset (there is an integer value from which the
number of records can be derived), or the datasets can be delimited by
a pattern (only the `name` field of the species and a leading '`-`' on
some numerical values occupies the first column). I prefer the latter
as species are then split out into distinct sets of records and can be
passed to a factory-type function/object to generate the internal
representation of the species.

So far we essentially have 3 broad categories of data, each containing
a number of species datasets. The next step is to decide how to
both parse the species datasets and represent them internally. This
necessitates a species dataset parser, and a class to facilitate
creation of these species. Finally, we need to have some sort of
framework for all of this to operate within. It makes sense to include
this as a dedicated module (the "get stuff from `thermo.inp`" module).
Within this module, a single point of access makes sense.

Bringing this all together leads to a set of different components to
the parser:

  - A file reader function to return the useful contents of the file,
    e.g. `read_thermoinp()`
  - A content-splitter function to return 3 distinct sets of data (gas
    products, condensed products, reactants), e.g. `_categorise()`
  - A class for representing chemical species, e.g. `Species` NOTE:
    this probably doesn't live in the same place as the rest of this.
  - A 'converter' for parsing per-species datasets and returning
	chemical species objects, e.g. `_create_species()`
  - This all lives in a module with a single point of access, e.g. the
	`thermoinp` module (named after the CEA-distributed data file)
	with a function `parse()` to return a dictionary of gaseous, 
	condensed and reactant species.

[CEA]: http://www.grc.nasa.gov/WWW/CEAWeb/index.htm
