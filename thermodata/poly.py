import collections

from thermodata import constants


# --------------------------------------------------------------------
# Classes
# --------------------------------------------------------------------
# Common elements of docstrings (cleared from namespace later)
_npdoc_fields = """
Fields
------

    lim : (T_min, T_max) [K]
    a : (a1, ..., an) coefficients
    b : (b1, b2) integration constants
    n : number of coefficients
    exp : exponent values
    dh : enthalpy diff.
"""

_npdoc_body = """
2002-spec with support for 8 terms, variable form and two constants
of integration. Contains polynomial coefficients for dimensionless
heat capacity at constant pressure Cp, where

    Cp/R' = p(T).

The polynomial may be integrated to provide dimensionless enthalpy,
H/R' and entropy S/R'T where the molar gas constant, R' is:

    R = {:.6f} J/(kmol K)
""".format(constants.R_CEA)


# Basic NASAPoly class; contains all metadata.
NASAPoly = collections.namedtuple('NASAPoly', 'lim, a, b, n, exp, dh')
NASAPoly.__doc__ = (
    """NASA Polynomial
    """ +
    _npdoc_body +
    _npdoc_fields
)


# NASAPoly with method to calculate non-dimensional heat capacity.
class NASAPolyND(NASAPoly):
    __doc__ = '\n'.join((
    "NASA Polynomial with non-dimensional quantity methods.",
    _npdoc_body,
    ("Methods are included to calculate the non-dimensional cp as a "
    "function of temperature."),
    _npdoc_fields
    ))

    def cpnd(self, T):
        """Return the non-dim. heat cap. at const. pressure [ND]."""
        return _dimless_heat_capacity(T, self.a)


# NASAPoly class with method to calculate molar heat capacity.
class NASAPolyML(NASAPoly):
    __doc__ = '\n'.join((
        "NASA Polynomial with molar quantity methods.",
        _npdoc_body,
        ("Methods are included to calculate the molar cp as a "
         "function of temperature."),
        _npdoc_fields
    ))

    def cpmol(self, T):
        """Return molar heat cap. at const. pressure [J/(kmol K)]."""
        return _dimless_heat_capacity(T, self.a) * constants.R_CEA


# --------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------
def _dimless_heat_capacity(T, a):
    # Returns the dimensionless heat capacity, Cp/R
    # T : Temperature, K
    # a : coefficients, len(a) == 7
    return (  a[0] / T**2
            + a[1] / T
            + a[2]
            + a[3] * T
            + a[4] * T**2
            + a[5] * T**3
            + a[6] * T**4
            )


# Tidy namespace
del _npdoc_body, _npdoc_fields
