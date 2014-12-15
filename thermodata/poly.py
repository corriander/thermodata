import math
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
        """Return non-dim. heat cap. at const. pressure []."""
        return _dimless_heat_capacity(T, self.a)

    def hnd(self, T):
        """Return non-dimensional enthalpy []."""
        return _dimless_enthalpy(T, self.a, b[0])

    def snd(self, T):
        """Return non-dimensional (T-dependent) entropy []."""
        return _dimless_entropy(T, self.a, b[1])


# NASAPoly class with method to calculate ndar heat capacity.
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

    def hmol(self, T):
        """Return molar enthalpy [J/kmol]."""
        return (_dimless_enthalpy(T, self.a, self.b[0]) *
                constants.R_CEA * T)

    def smol(self, T):
        """Return molar (T-dependent) entropy [J/(kmol K)]."""
        return _dimless_entropy(T, self.a, self.b[1]) * constants.R_CEA


# Parser class for handling datasets -> (NASAPoly*, ...)
class Parser(object):

    def __init__(self, polycls=NASAPoly):
        self.polycls = polycls

    def __call__(self, dataset):
        """Return NASAPoly instances for each interval in a dataset.
        """
        # Split dataset into records
        records = dataset.strip().split('\n')

        # No. intervals is at (2, 2) or (1, 1) in 0-index
        if records[1][1] == 0:
            polys = ()
        else:
            polys = self._parse_intervals(records[2:])
        return polys

    def _parse_intervals(self, records):
        # Return a tuple of NASAPoly* instances for a list of records
        # containing interval metadata and polynomial specification.

        return tuple(
            self._parse_interval(lines)
            for lines in self._intervals(records)
        )

    @staticmethod
    def _intervals(records):
        # Return records split into triplets
        for i in range(0, len(records), 3):
            yield records[i:i+3]

    def _parse_interval(self, records):
        # Return a NASAPoly* instance for a record triplet

        # the first line is metadata, the second two specify the poly
        lim, n, exp, dh = self._parse_metadata(records[0])
        a, b = self._parse_coefficients(records[1:])

        return self.polycls(lim, a, b, n, exp, dh)

    def _parse_metadata(self, record):
        # Returns a tuple (lim, n, exp, dh)
        lim = tuple(float(n) for n in record[:22].split())
        n = int(record[22])
        exp = tuple(float(n) for n in record[23:63].split())
        dh = float(record[65:])

        return (lim, n, exp, dh)

    def _parse_coefficients(self, records):
        # Returns tuple of coefficients (a1, .. a7) & consts (b1, b2)
        array1, array2 = records

        # parse records containing numerical strings
        a = self._double_array_to_float(array1)
        a.extend(self._double_array_to_float(array2[:32]))
        a = tuple(a)
        b = tuple(self._double_array_to_float(array2[48:]))

        return (a, b)

    def _double_array_to_float(self, string):
        # Parse a string a containing 16-char Fortran-style doubles into
        # a list of floats
        float_strings = [string[i:i+16].replace('D','e') # Pythonify
                         for i in range(0, len(string), 16)]
        return list(map(float, float_strings))


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

def _dimless_enthalpy(T, a, b):
    # Returns the dimensionless enthalpy, H/RT
    # T : Temperature, K
    # a : coefficients, len(a) == 7
    # b : integration constant
    return (- a[0] / T**2
            + (a[1] * math.log(T) + b) / T
            + a[2]
            + a[3] * T / 2.0
            + a[4] * T**2 / 3.0
            + a[5] * T**3 / 4.0
            + a[6] * T**4 / 5.0
            )

def _dimless_entropy(T, a, b):
    # Returns the dimensionless entropy, S/R.
    # T : Temperature, K
    # a : coefficients, len(a) == 7
    # b : integration constant
    return (- a[0] / T**2 / 2.0
            - a[1] / T
            + (a[2] * math.log(T) + b)
            + a[3] * T
            + a[4] * T**2 / 2.0
            + a[5] * T**3 / 3.0
            + a[6] * T**4 / 4.0
            )


# Tidy namespace
del _npdoc_body, _npdoc_fields
