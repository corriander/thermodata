"""Interface for NASA Glenn thermodynamic database source file.

This module provides low-level access to the original database format.
The database is broadly split into three categories of chemical
species; gaseous equilibrium products, condensed equilibrium products
and reactants (e.g. 'Air' or other mixtures with characteristic
properties).

The DB class is the main point of access now. It provides categories
of species data (as SpeciesRecords) amongst other functionality.
"""
import re
import os
import collections

from thermodata import poly


class DB(object):
    """Interface to the 'thermo.inp' source database.

    The source database is grouped into three distinct sections:

        condensed products and reactants
        gaseous products and reactants
        reactants

    This class provides subsets of the database species.
    """

    polytype = poly.NASAPoly

    # Define the thermo.inp format header.
    header = '\n'.join([
        '{:<80s}'.format('thermo'),
        '   200.000  1000.000  6000.000 20000.000   9/09/04'
    ])

    def __init__(self, polytype=''):
        self._select_polytype(polytype)
        self._parse()
        self._dict = {s.name:s for s in self.all}

    # ----------------------------------------------------------------
    # Categories
    # ----------------------------------------------------------------
    # TODO: Categories now unimportant with the functionality of
    # subset. They should be removed.
    @property
    def all(self):
        """All species."""
        return self._condensed + self._gaseous + self._reactant

    @property
    def allcondensed(self):
        """All condensed species (including reactants)."""
        return [s
                for s in (self._condensed + self._reactant)
                if s.phase > 0]

    @property
    def allgases(self):
        """All gaseous species (including reactants)."""
        return [s
                for s in (self._gaseous + self._reactant)
                if s.phase == 0]

    @property
    def product(self):
        """Species that can appear as products in reactions."""
        return self._condensed + self._gaseous

    # NOTE wrt above TODO: Keep these!
    @property
    def condensed(self):
        """Condensed, product-only species."""
        return self._condensed

    @property
    def gaseous(self):
        """Gaseous, product-only species."""
        return self._gaseous

    @property
    def reactant(self):
        """Mixed-phase, reactant-only species."""
        return self._reactant

    # ----------------------------------------------------------------
    # External methods
    # ----------------------------------------------------------------
    def format(self):
        """Return database in syntactically valid string format.

        The string should be parsable by applications able to read
        the source format.

        Example
        -------

            >>> db = DB() # empty database
            >>> print(db.formatted)
            thermo
               200.000  1000.000  6000.000 20000.000   9/09/04
            END PRODUCTS
            END REACTANTS
        """
        db = [self.header]
        db.append('\n'.join(s.formatted for s in self.condensed))
        db.append('\n'.join(s.formatted for s in self.gaseous))
        db.append('{:<80s}'.format('END PRODUCTS'))
        db.append('\n'.join(s.formatted for s in self.reactant))
        db.append('{:<80s}'.format('END REACTANTS'))

        return '\n'.join(filter(None, db))

    def list_categories(self):
        """List categories implemented in the original database."""
        return ['condensed', 'gaseous', 'reactant']

    def list_species(self, category=''):
        """List species in the database.

        Arguments
        ---------

            category : filter by category
        """
        l = []
        if category:
            categories = [category,]
        else:
            categories = self.list_categories()

        for category in categories:
            l.extend([s.name for s in getattr(self, category)])

        return l

    def lookup(self, string):
        """Query the database for species names matching `string`.

        Returns a list of SpeciesRecord.

        Usage
        -----

        This method wraps re.match so regexen works.

            >>> lst = db.lookup('*.H2')
            >>> len(lst)
            78

        Note that parantheses are present in many species names and
        are therefore escaped by default. Don't treat parentheses in
        any special way.

            >>> len(db.lookup('Jet-A(g)'))
            1
            >>> len(db.lookup('Jet-A\(g\)'))
            0

        Batteries aren't included for filtering by category, but it's
        simple enough to construct more advanced queries. For example,
        searching for all species with 'H2' in the name that are also
        reactants might look a little like:

            >>> [s.name for s in db.lookup('.*H2') if s in db.reactant]
            ['(CH2)x(cr)', 'C2H2(L),acetyle', 'C6H5NH2(L)', 'H2(L)', 'H2O2(L)']
        """
        string = string.replace('(','\(').replace(')','\)')
        lst = []
        for name, obj in self._dict.items():
            if re.match(string, name):
                lst.append(obj)
        return sorted(lst, key=lambda o: o.name)

    def subset(self, species=(), filt=None):
        """Create a subset of this database.

        Returns a new DB instance containing species matching the
        criteria. There are two methods of defining criteria, by a
        species name pattern (or iterable of patterns) which get
        passed to `lookup` and a filter function.

        Arguments
        ---------

            species : string or iterable of strings containing
                patterns which get matched against species names.
            filt : callable that takes an object (SpeciesRecord) and
                returns a boolean value. Used directly in `filter`.

        Examples
        --------

        Filter function for any gas species:

            >>> subset = DB().subset(filt=lambda o: o.phase == 0)

        Species list:

            >>> subset = DB().subset(('.*H2', 'Air'))
        """
        # Argument handling
        # -----------------
        # We should accept a string here (e.g. '.*H2')
        if isinstance(species, str):
            species = (species,)

        # Create the None filter
        if filt is None:
            filt = lambda obj: True

        # Create a set of species matching the species specification
        if species:
            species_set = set()
            for string in species:
                species_set.update(set(self.lookup(string)))
        else:
            species_set = set(self._dict.values())

        # New DB instance (super crude, might want to subclass)
        # -----------------------------------------------------
        subset = self.__class__()
        subset._dict = {}
        for obj in filter(filt, species_set):
            subset._dict[obj.name] = obj

        # Bit of a hack workaround; previously the database was split
        # into three groups as it was parsed. This information is now
        # carried with the species (isproduct and phase). This is
        # fairly deeply entrenched in the database structure currently
        # so we're going to overwrite the condensed, gaseous and
        # reactant lists which are the guts (but _dict should be now)
        objs = subset._dict.values()
        reactants = filter(lambda o: o in self.reactant, objs)
        gaseous = filter(lambda o: o in self.gaseous, objs)
        condensed = filter(lambda o: o in self.condensed, objs)
        sort_key = lambda o: o.name
        subset._reactant = sorted(reactants, key=sort_key)
        subset._condensed = sorted(condensed, key=sort_key)
        subset._gaseous = sorted(gaseous, key=sort_key)

        return subset

    # ----------------------------------------------------------------
    # Internal methods
    # ----------------------------------------------------------------
    def _parse_to_categories(self):
        """Split database file into categories.

        File location is hardcoded.
        """
        categ_dict = _read_categories()
        self._condensed = categ_dict['condensed_products']
        self._gaseous = categ_dict['gas_products']
        self._reactant = categ_dict['reactants']

    def _parse_category(self, category):
        """Split category into species datasets."""
        pattern = re.compile(r'\n(?=[eA-Z(])')
        name = '_{}'.format(category)

        # Add a flag to indicate whether the species is reactant-only
        if category == 'reactant':
            isproduct = False
        else:
            isproduct = True

        # Split category (string) into species dataset (strings) and
        # cast them as SpeciesRecord instances.
        # FIXME: src should be passed directly, but for now other
        # functions in this module are dependent on this list form.
        l = []
        for src in pattern.split(getattr(self, name)):
            src = src.split('\n')
            polycls = self.polytype
            sr = SpeciesRecord.from_dataset(src, isproduct, polycls)
            l.append(sr)
        setattr(self, name, l)

    def _parse(self):
        """Split database file into (categorised) datasets."""
        self._parse_to_categories()

        for c in self.list_categories():
            self._parse_category(c)

    def _select_polytype(self, polytype):
        # Selects appropriate class from module: poly
        cls = getattr(poly, 'NASAPoly{}'.format(polytype.upper()))
        self.polytype = cls

    # ----------------------------------------------------------------
    # Magic methods
    # ----------------------------------------------------------------
    def __getitem__(self, key):
        """Retrieve SpeciesRecord by species name (dict-like)."""
        # Access the hidden _dict which is keyed with species name.
        return self._dict[key]


# --------------------------------------------------------------------
#
# Internal functions
#
# --------------------------------------------------------------------
def _pprint_refcode(code):
    """Look up the NASA GRC reference code and format."""

    references = {
        'g' : 'Glenn Research Center',
        'j' : 'NIST-JANAF Thermochemical Tables. Chase,1998',
        't' : ('Thermodynamic Properties of Individual Substances. '
               'Gurvich 1978, 1979, 1982, 1989, 1991, 1996'),
        'n' : 'TRC Thermodynamic Tables, NIST',
        'b' : 'Thermochemical Data of Pure Substances. Barin 1989',
        'c' : 'CODATA Key Values for Thermodynamics. Cox 1989',
        's' : 'Standard Reference Data: J.Phys.Chem.Ref.Data'
        }
    months = ('', 'Jan', 'Feb', 'Mar', 'April', 'May', 'Jun',
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
    expandyear = {'0' : '20'}
    regex = re.compile(r'\d.*')

    # Evaluate references
    reference = references.get(code[0])
    if code == 'g tpis':
        # Special case for C4
        reference += '{!s}, {!s}'.format(reference,
                                         references.get('t'))
        return '{!s}, 0000'.format(reference)

    # Evaluate date
    date_code = regex.search(code).group()
    date = []
    try:
        month, year = date_code.split('/')
        date.extend((expandmonth[int(month)], '. '))
    except ValueError:
        year = date_code
    date.extend((expandyear.get(year[0], '19'), year))

    return 'Reference       : {!s}\nDate Calculated : {!s}'.format(
            reference,
            ''.join(date)
            )

# TODO: To be deprecated - replace with soft-coded file approach
def _read_categories():
    # Split the database into three category strings.
    # Returns a category-keyed dictionary of string values.
    path = os.path.join(os.path.dirname(__file__),
                        'data',
                        'thermo.inp')
    keys = 'gas_products', 'condensed_products', 'reactants'
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
    return dict(zip(keys, pattern.search(contents).groups()))

# --------------------------------------------------------------------
#
# Data structure for the 2002-spec polynomial form (variable form
# <=8-term poly with 2 integration constants)
#
# --------------------------------------------------------------------
_Interval = collections.namedtuple('TemperatureInterval',
                                   ['bounds',
                                    'ncoeff',
                                    'exponents',
                                    'deltah',
                                    'coeff',
                                    'const'])



class Interval(_Interval):
    """Specification of polynomial function of temperature.

    This class of object stores data describing a variable form
    polynomial function applicable to a defined temperature interval.
    Fields correspond to the following data:

      `bounds`    : Interval bounds (Tmin, Tmax)
      `ncoeff`    : Number of coefficients/terms
      `exponents` : Exponent magnitudes (len() == ncoeff)
      `deltah`    : Reference enthalpy value
      `coeff`     : Coefficients (len() == ncoeff)
      `const`     : Integration constants

    """
    pass

def _double_array_to_float(string):
    # Parse a string a containing 16-char Fortran-style doubles into
    # a list of floats
    float_strings = [string[i:i+16].replace('D','e') # Pythonify
                     for i in range(0, len(string), 16)]
    return list(map(float, float_strings))

def _parse_interval(records, cls=Interval):
    # Parse records containing a temperature interval/polynomial spec.
    # This expects records as a list of strings and returns an
    # Interval instance.
    metadata, array1, array2 = records

    # parse metadata string first
    bounds = tuple(float(n) for n in metadata[:22].split())
    ncoeffs = int(metadata[22])
    exponents = tuple(float(n) for n in metadata[23:63].split())
    deltah = float(metadata[65:])

    # parse records containing numerical strings
    coeffs = _double_array_to_float(array1)
    coeffs.extend(_double_array_to_float(array2[:32]))
    coeffs = tuple(coeffs)
    consts = tuple(_double_array_to_float(array2[48:]))

    if issubclass(cls, poly.NASAPoly):
        args = bounds, coeffs, consts, ncoeffs, exponents, deltah
    else:
        args = bounds, ncoeffs, exponents, deltah, coeffs, consts

    return cls(*args)

# --------------------------------------------------------------------
#
# Data structure for the 2002-spec species datasets
#
# --------------------------------------------------------------------
_Species = collections.namedtuple('SpeciesRecord',
                                  ['name',
                                   'comments',
                                   'nintervals',
                                   'refcode',
                                   'formula',
                                   'phase',
                                   'molwt',
                                   'h_formation',
                                   'h_assigned',
                                   'T_reference',
                                   'intervals',])

# TODO: Make intervals hashable, currently is a list.
class SpeciesRecord(_Species):
    """Chemical species metadata and thermodynamic properties.

        `name` 		  : Species name/ID (usually formula)
        `comments`    : References and comments
        `nintervals`  : Number of temperature intervals (0-3)
        `refcode`     : Reference-date code
        `formula`     : String array of elements and no. atoms
        `phase`       : Condensed phases non-zero
        `molwt`       : molecular weight/molar mass, kg/kmol
        `h_formation` : heat/enthalpy of formation (nintervals > 0)
        `h_assigned`  : assigned enthalpy (nintervals == 0)
        `T_reference` : reference temperature (for assigned enthalpy)
        `intervals`   : temperature intervals (nintervals > 0)
    """

    @property
    def formatted(self):
        """Return species dataset as a thermo.inp formatted string."""
        return self._formatted

    @property
    def isproduct(self):
        """Flag indicates if species is a valid reaction product."""
        return self._isproduct


    @classmethod
    def from_dataset(cls, records, isproduct=False, polycls=Interval):
        """Create a SpeciesRecord instance from a thermo.inp block.

        Arguments
        ---------

            records : list of thermo.inp species records (strings)

        Each species dataset has a number of records/lines (3-11)
        """
        # Parse records containing species data.
        # Returns a Species instance.

        # We want to keep the source data around
        # FIXME: this undoes a previous operation.
        rawdata = '\n'.join(records)

        # split the records up
        head, body, tail = records[0], records[1], records[2:]

        # Parse the name & comments from the header
        name, comments = _parse_first_record(head)

        # Parse the non-polynomial data
        nintervals = int(body[1])
        refcode = body[2:10].strip()
        # make formula a bit more parse-friendly but leave as a string
        # e.g.
        #	'C   1.00O  2.00   0.00   0.00   0.00' -> 'C:1.00 O:2.00'
        formula = ' '.join([
            '{!s}:{!s}'.format(body[i:i+2].strip(), body[i+2:i+8].strip())
            for i in range(10, 50, 8)
            ]).replace(' :0.00', '')
        phase = int(body[51])
        molwt = float(body[52:65])

        # At this stage, the fields have potential to vary depending
        # on whether temperature intervals are present or not.
        refenthalpy = float(body[65:])
        if nintervals > 0:
            h_assigned = T_reference = None
            h_formation = refenthalpy
            # each interval is described by three records
            intervals = tuple(_parse_interval(tail[i:i+3], polycls)
                              for i in range(0, len(tail), 3))
        else:
            # FIXME: intervals should probably be an empty tuple.
            h_formation = intervals = None
            h_assigned = refenthalpy
            T_reference = float(tail[0].split()[0]) # grab first word

        inst = cls(name,
                   comments,
                   nintervals,
                   refcode,
                   formula,
                   phase,
                   molwt,
                   h_formation,
                   h_assigned,
                   T_reference,
                   intervals)

        inst._formatted = rawdata
        inst._isproduct = isproduct
        return inst


def _parse_species(records):
    return SpeciesRecord.from_dataset(records)

def _parse_first_record(record):
    # Takes the first record of a species dataset and returns the name
    # and comment fields
    return record[:18].rstrip(), record[18:].rstrip()

