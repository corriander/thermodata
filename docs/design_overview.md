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

[CEA]: http://www.grc.nasa.gov/WWW/CEAWeb/index.htm
