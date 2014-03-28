Design Overview
===============

This code is intended primarily to provide an alternative,
Python-based interface to the NASA Glenn thermodynamic database
distributed with the NASA Glenn program [Chemical Equilibrium with
Applications (CEA)][CEA] by Sanford Gordon and Bonnie J. McBride. A
secondary goal is to wrap CEA itself. This document summarises the
rationale and high-level design approach to achieving this. For
details about the database, its purpose and usage, see <database.md>.

Rationale
---------

This code is a byproduct of my research a few years ago. Originally I
was using CEA to perform calculate equilibrium chemical composition
but required access to a small subset of the database to perform
external calculations. [ThermoBuild][] provides web-based support for
generating both tabular data *and* subsets of the database. However,
it's not suitable for automated access, especially for variable
subsets. I wrote a simple `grep`-based utility (if you can go as far
as calling it that) for pulling arbitrary species datasets but I
always intended to provide a more general method of getting at the
data (including parsing it) but never had the time to implement it.
This is, in part, a realisation of that intent.

Regarding wrapping CEA: Accessing the program in an automated manner
is possible via your chosen system shell. This is fine for simple use
cases and even scales up to repeated execution. I'd already
implemented this in a couple of different ways in Scilab, so this is
basically a port at the simplest level but doesn't suffer from the
same language limitations.

Design Specification
--------------------

The following is a loose design specification oriented around the two
goals of interfacing with the database and wrapping CEA.

### Database Interface

There are several components to this:

  - A source database parser.
  - An alternative database format (e.g. XML) including facilities for
    writing and parsing/querying.
  - Searching/browsing functionality.
  - Methods for creating predefined and/or custom subsets of the data.

The database parser should support the full database
specification including redundant information. The alternative format
need not support this information, though it should be easily
extendable.

Browsing/searching involves creating suitable species representations
and methods for searching species name and comment data. Although
these functions are somewhat independent of the backend, it may be
helpful to support retrieval of records from the source format as well
as in native Python data structures (a la [ThermoBuild][]).

### CEA Interface

Simply, the CEA interface needs to provide a mechanism for preparing
input files, reading output files and calling the CEA program
efficiently and cleanly. This is relatively trivial, although could be
made more sophisticated by circumventing the need for disk I/O.

[CEA]: http://www.grc.nasa.gov/WWW/CEAWeb/index.htm
