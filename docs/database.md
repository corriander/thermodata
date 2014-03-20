CEA Database Overview
=====================

The CEA thermo database used here has data for over 2000 chemical
species. Species are broadly grouped into two types, products and
reactants (the latter only existing as a reactant not as an
equilibrium product). However, within these categories species may be
in a solid, liquid or gaseous state. Where a chemical species is
present in different states or in ionic forms, these are treated as
distinct species (e.g. the element hydrogen and its positive and
negative ions is represented as 3 distinct species: `H`, `H+` and
`H-`).

The species data is in the form of empirical coefficients for a
polynomial function of temperature (capable of calculating specific
heat capacity, enthalpy and entropy at a given temperature). These
coefficients are provided for up to 3 temperature intervals (by and
large, 200 - 20,000 K) though for some reactant species no interval is
provided.

The database format is ASCII, with 80-char width records split into
fields of defined (field-dependent) widths. 

Fields
------

For a given species there are between 2 and 11 records depending on
the number of provided temperature intervals (0-3).

### Record 1

  -	Species name/formula (chars 1-15). This serves as an ID. Note that
  	'l' is represented by `L` and condensed phases designated as α, β,
	γ or δ are renamed *a*, *b*, *c* or *d* due to ASCII limitations.
  - Comments (chars 16-80). These include references in the format of
    author, year or page and date in the case of TRC tables. When heat
	of formation is taken from a separate reference, this is included 
	as `Hf:<ref>`. Reference elements or species used for heat of
	formation calculations are indicated by `Ref-Elm` or
	`Ref-Species`.

### Record 2

  - Reference-Date code (1-6)

### Record 3

  - Temperature interval bounds min/max 

Redundancy
----------

The current version of the database has built-in redundancy. 

  - The polynomials use 7 coefficients but there is room for 8. 
  - Variable temperature intervals are supported (and used for
    non-gaseous species).
  - The polynomial can be defined in terms of variable, non-integer
    temperature exponents. In practice, the exponents are both integer
	and common across all species.

States
------

All data is provided for species in their standard state at the
specified temperature. For gases, this is the standard-state pressure
of 1 bar (100,000 Pa). For condensed species this is the standard
atmosphere pressure, 1 atm (101,325 Pa).

For chemical elements, the reference state is the stable state at
298.15 K. For gaseous elements at this temperature (and the
standard-state pressure), the reference state is a gas over the entire
temperature range. Condensed elements are treated as condensed over
the temperature range with phase-transitions. These are identified by
"Ref-elm" or "Ref-species" in the comments.

The universal gas constant is defined as *R* = 8.314510 J/mol/K.

References
----------

  - Gordon, McBride and Zehe *NASA Glenn Coefficients for Calculating
    Thermodynamic Properties of Individual Species*. Technical Report
	NASA/TP-2002-211556 2002.
