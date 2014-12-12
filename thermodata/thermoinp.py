"""Interface for NASA Glenn thermodynamic database source file.

This module provides low-level access to the original database format.
The database is broadly split into three categories of chemical
species; gaseous equilibrium products, condensed equilibrium products
and reactants (e.g. 'Air' or other mixtures with characteristic
properties).

The DB class is the main point of access now. It provides categories
of species data (as SpeciesRecords) amongst other functionality.

Deprecated:

Several API functions are included to extract data from the source
file at several levels of decomposition.

  - `parse` returns a category-keyed dictionary of species datasets
    parsed into a namedtuple.
  - `lookup` can search the database for species with a name matching
    a prefix.
  - `list_species` provides simple categorised lists of species names.

"""
import re
import os
import collections

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

# NOTE: Deprecated (DB._parse_category)
def _read_species():
    # Split the database into categorised lists of species.

    # Returns a category-keyed dictionary of species-separated
    # content.  Each category is a list of species dataset strings.
    # The strings contain newline-delimited records containing
    # species-specific data.
    categories = _read_categories()
    # For each category, separate into per-species strings
    pattern = re.compile(r'\n(?=[eA-Z(])')
    return {k:pattern.split(v) for k, v in categories.items()}

# NOTE: Deprecated (DB._parse)
def parse():
    """Parse the database into categorised lists of Species instances.

    Returns a category-keyed dictionary of lists containing Species
    instances.

        >>> d = parse()
        >>> sorted(d.keys())
        ['condensed_products', 'gas_products', 'reactants']
        >>> type(d['reactants'])
        <type 'list'>
        >>> sample_species = d['reactants'][0]
        >>> sample_species.name
        'Air'

    """
    species_categories = _read_species()
    return {k:[_parse_species(string.split('\n')) for string in lst]
            for k, lst in species_categories.items()
            }

# NOTE: Deprecated (DB.__getitem__)
def lookup(prefix, form='parsed', exact=False):
    """Locate species with a matching name prefix.

    This will search the database and return species datasets where
    the name matches a prefix, returning the results in a
    category-keyed dictionary of lists (no results in a category is
    represented by None value).

        >>> matches = lookup('N2')
        >>> type(matches)
        <type 'dict'>
        >>> len(matches['reactants'])
        3
        >>> matches['reactants'][1].name
        'N2H4(L)'

    Optionally, the search can return unprocessed species datasets:

        >>> matches = lookup('N2', form='unparsed')
        >>> matches['reactants'][1][:40] # limit sample string length
        'N2H4(L)           Hydrazine.Hf:Gurvich,1'

    The search can be category-restricted neatly using Python syntax:

        >>> matches = lookup('N2')['reactants']
        >>> type(matches)
        <type 'list'>
        >>> len(matches)
        3

    More advanced searching/browsing is outside the scope of this
    module.

    """
    if exact:
        pattern = re.compile(r'{}$'.format(re.escape(prefix)))
    else:
        pattern = re.compile(r'{}'.format(re.escape(prefix)))

    if form == 'parsed':
        # parse the database file into Species objects, set up the
        # match function to check against the species name.
        source = parse()
        def match(species):
            return pattern.match(species.name)

    elif form == 'unparsed':
        # split the database file into per-species strings, set up the
        # match function to check against the name region of the
        # string (with trailing whitespace stripped)
        source = _read_species()
        def match(species):
            name, __ = _parse_first_record(species.split('\n')[0])
            return pattern.match(name)

    else:
        raise ValueError("Argument form='{!s}' invalid".format(form))

    results = dict.fromkeys(source)	# results container
    for category in results:
        # keep track of matches
        matches = [species
                   for species in source[category]
                   if match(species)
                   ]
        if matches:
            # if there are matches, add them to results
            results[category] = matches

    return results


class DB(object):
    """Interface to the 'thermo.inp' source database."""

    # Define the thermo.inp format header.
    header = '\n'.join([
        '{:<80s}'.format('thermo'),
        '   200.000  1000.000  6000.000 20000.000   9/09/04'
    ])

    def __init__(self):
        self._parse()
        self._dict = {s.name:s for s in self.all}

    # ----------------------------------------------------------------
    # Categories
    # ----------------------------------------------------------------
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
        db.append('END PRODUCTS')
        db.append('\n'.join(s.formatted for s in self.reactant))
        db.append('END REACTANTS')

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
        l = [SpeciesRecord.from_dataset(src.split('\n'), isproduct)
             for src in pattern.split(getattr(self, name))]
        setattr(self, name, l)

    def _parse(self):
        """Split database file into (categorised) datasets."""
        self._parse_to_categories()

        for c in self.list_categories():
            self._parse_category(c)

    # ----------------------------------------------------------------
    # Magic methods
    # ----------------------------------------------------------------
    def __getitem__(self, key):
        """Retrieve SpeciesRecord by species name (dict-like)."""
        # Access the hidden _dict which is keyed with species name.
        return self._dict[key]

# TODO: To be deprecated: represent categories by set in DB
def create_subset(prefix=None,
                  filter_category=None,
                  exact=False):
    """Return a syntactically valid subset of thermo.inp in a string.

    This function provides means to generate a custom database in the
    original format that should be readable by applications that parse
    it directly (e.g. CEA). See below for an explanation of the
    general specification for this format.


    Arguments
    ---------

        prefixes : search string prefix (uses re.match())
        filter_category : restrict search to category:

             'gas_products',
             'condensed_products',
             'reactants'

        exact : performs a whole-word search rather than prefix match


    Usage
    -----

    The prefix is passed to lookup() and exhibits the same behaviour.
    For example, to create a subset containing the family of jet
    fuels:

        >>> string = create_subset('JP')

    To create a subset of all gaseous products with names containing
    the prefix 'N2':

        >>> string = create_subset('N2', 'gas_products')

    To output all reactant species:

        >>> string = create_subset(filter_category='reactants')

    Note that the lookup is based on matching a species name with a
    prefix. The exact parameter can be specified to match whole name
    strings.


    Database Format
    ---------------

    The general structure of the thermo.inp data format is:

        thermo
           200.000  1000.000  6000.000 20000.000   9/09/04
        # Product species datasets
        END PRODUCTS
        # Reactant species datasets
        END REACTANTS

    As long as these features are present the database is
    syntactically valid and should be parsable.
    """

    # Recreate the delimiters, these can be interleaved with
    # categories. This is quicker than pulling them from the source
    # and easier than passing this data around.
    header = ['{:<80s}'.format('thermo')]
    intervals = '   200.000  1000.000  6000.000 20000.000   9/09/04'
    header.append(intervals)
    delimiters = ('\n'.join(header),
                  '',
                  '{:<80s}'.format('END PRODUCTS'),
                  '{:<80s}'.format('END REACTANTS')
                  )

    # empty list for blocks of species, map category to an index
    matching_species = [[], [], []] # 3 categories of species
    index = {'gas_products' : 0,
             'condensed_products' : 1,
             'reactants' : 2
             }

    if (prefix, filter_category) == (None, None):
        # this is a no-op
        raise ValueError("No subset criteria selected")

    elif prefix is None:
        # i.e. filter_category is not None
        # category of species can be obtained directly as a string
        string = _read_categories()[filter_category]
        matching_species[index[filter_category]] = string

    else:
        # prefix is defined. category might be.
        if isinstance(prefix, str):
            prefix = [prefix]
        for string in prefix:
            category_dict = lookup(string, 'unparsed', exact)

            if filter_category is not None:
                # We can filter the lookup results, giving us a list.
                # Reference the correct element of the
                # matching_species container
                dataset_list = category_dict[filter_category]
                output_list = matching_species[index[filter_category]]
                if dataset_list is not None:
                    # Extend the list with new datasets.
                    output_list.extend([dataset
                                        for dataset in dataset_list
                                        if dataset not in output_list]
                                       )
                continue

            # ^^ continue means we're only here if category is not
            # specified
            for category, dataset_list in category_dict.items():
                # Start by referencing the appropriate element of the
                # matching_species container
                output_list = matching_species[index[category]]
                if dataset_list is not None:
                    # Extend the list with new datasets
                    output_list.extend([dataset
                                        for dataset in dataset_list
                                        if dataset not in output_list]
                                       )

        # All search strings have been looked up, matching_species is
        # a list of three lists containing matching datasets. These
        # need processing.
        for i, list_ in enumerate(matching_species):
            sorted_datasets = sorted(list_)
            matching_species[i] = '\n'.join(sorted_datasets)

    subset = [None] * 7 # len(delimiters) + len(species)
    subset[::2] = delimiters # populate even with delimiters
    subset[1::2] = matching_species # fill odd indices with matches

    return '\n'.join(filter(None, subset))

# NOTE: Deprecated
def list_species():
    """List species in the database.

    Returns a category-keyed dictionary of species names

    """
    return {category:[species.name for species in species_list]
            for category, species_list in parse().items()}

# --------------------------------------------------------------------
#
# General utility functions (non-API)
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
    def from_dataset(cls, records, isproduct=False):
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
            intervals = [_parse_interval(tail[i:i+3])
                         for i in range(0, len(tail), 3)]
        else:
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

def _parse_interval(records):
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

    return Interval(bounds,
                    ncoeffs,
                    exponents,
                    deltah,
                    coeffs,
                    consts)
