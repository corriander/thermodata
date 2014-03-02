python-thermo
=============

Thermodynamic data, functions and objects in Python

This repository is a collation/re-working of data structures,
thermodynamic data and functions originally written in Scilab and
ported to Python (with the added benefit of some hindsight).

Thermodynamic Data
------------------

Currently, this package provides a Python representation of the
thermodynamic dataset for over 2000 chemical species distributed with
the venerable [Chemical Equilibrium with Applications][CEAWeb] FORTRAN
code developed by Bonnie J. McBride and Sanford Gordon. The `cea`
package will probably provide a wrapper for this code at some point,
but for now it's just an alternative database access point.

Issues
------

  -	I really would like some sort of tests for the integrity of the
  	`thermo.inp` data after being pulled into object representation.
  - Reactants and Products are currently not distinguished in this
    dataset.

[CEAWeb]: http://www.grc.nasa.gov/WWW/CEAWeb/ceaHome.htm

