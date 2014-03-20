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

The database format is ASCII, with 80-column records split into
fields of defined (field-dependent) widths. 

Fields
------

For a given species there are between 2 and 11 records depending on
the number of provided temperature intervals (0-3). Each temperature
interval has 3 records, detailed here in Record 3, 4 and 5. For
multiple temperature intervals (nonionic gases have 2, or 3 in the
case of monatomic species and simple molecules) these records are
repeated for each interval.

Included in this is shorthand for the field columns and data type. In
the case of strings and integers it's fairly obvious (int is a 1 digit
integer, 2int is a 2-digit integer, 15str is a 15-digit string, etc).
For floats, the shorthand mirrors pythonic formatting, so a 10.3f is a
float comprised of a 10-digit integer and a 3-digit decimal component.
Note that in the database, FORTRAN 77 differentiates between floats
and doubles whereas python doesn't. Therefore, doubles have been
described here as floats.

### Record 1

  -	**Species name/formula (cols 1-15, 15str)**. This serves as
  	an ID. Note that 'l' is represented by `L` and condensed phases
	designated as α, β, γ or δ are renamed *a*, *b*, *c* or *d* due to
	ASCII limitations.
  - **Comments (cols 16-80, 65str)**. These include references in
    the format of author, year or page and date in the case of TRC
	tables. When heat of formation is taken from a separate reference,
	this is included as `Hf:<ref>`. Reference elements or species used
	for	heat of formation calculations are indicated by `Ref-Elm` or
	`Ref-Species`.

### Record 2

  - **Number of temperature intervals (col 2, 2int)**.
  - **Reference-Date code (cols 4-9, 6str)**. This includes a
    character indicating a general reference followed by a date (e.g.
	`g` indicates that NASA Glenn was the source of significant work 
	in deriving the data and `10/96` indicates the month/year).
  - **Chemical formula (cols 11-50, 2str + 6.2f)**. This is a set of 5 
	element/atom, number pairs. In the vast majority of cases the 
	numbers are integers but in some cases they are non-integer, so 
	floats are used.
  - **Phase (col 52, int)**. Zero for gas, nonzero for condensed phases.
  - **Molar mass (cols 53-65, 13.5f)**. Originally labelled molecular
    weight (in units g/mol).
  - **Heat of formation (cols 66-80, 13.5f)**. In the case of condensed
    species this is actually an assigned enthalpy (equivalent to the
	heat of formation at 298.15 K). Units J/mol.

### Record 3

  - **Temperature range (cols 2-21, 2x 10.3f)**. The minimum and maximum
    bounds for the current temperature interval. Units, K.
  - **Number of coefficients (col 23, int)**. This is always 7 in this
    data (though the database format supports 8, see section
	Redundancy).
  - **Polynomial exponents (cols 24-63, 8x 5.1f)**. These are always `[-2,
    -1, 0, 1, 2, 3, 4]` in this data.
  - **{H(298.15) - H(0)} (cols 66-80, 15.3f)**. This is the difference
    between the heat of formation at the enthalpy at T = 0 K.

### Record 4

  - **Coefficients 1-5 (cols 1-80, 5x 16.8f)**.

### Record 5

  - **Coefficients 6-8 (cols 1-48, 3x 16.8f)**. The 8th is not used in
    this data (see section Redundancy).
  - **Integration constants (cols 49-80, 2x 16.8f)**. Used in evaluation
	of enthalpy and temperature-dependent component of entropy,
	respectively.


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
`Ref-elm` or `Ref-species` in the comments.

The universal gas constant is defined as *R* = 8.314510 J/mol/K.

References
----------

  - Gordon, McBride and Zehe *NASA Glenn Coefficients for Calculating
    Thermodynamic Properties of Individual Species*. Technical Report
	NASA/TP-2002-211556. 2002.
  - Gordon and McBride *Computer Program for Calculation of Complex
    Chemical Equilibrium Compositions and Applications Part II*.
	NASA-RP1311b. 1996.
