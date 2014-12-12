import os
import unittest
import thermodata.thermoinp as thermoinp

Species = thermoinp.SpeciesRecord
_Interval = thermoinp._Interval


class TestParse(unittest.TestCase):
    def setUp(self):
        self.parsed = thermoinp.parse()
    
    def test_gas(self):
        """Test sample gas in results.

        The test data has been constructed independently from the code
        data flow, a match here implies that the data is being
        parsed accurately.
        
        """
        self.assertTrue(test_gas in self.parsed['gas_products'])
    
    def test_condensed(self):
        """Test sample condensed products in results.

        The test data has been constructed independently from the code
        data flow, a match here implies that the data is being
        parsed accurately.
        
        """
        self.assertTrue(test_condensed in 
                        self.parsed['condensed_products']
                        )
    
    def test_reactant(self):
        """Test sample reactants are in the reactants category.


        The test data has been constructed independently from the code
        data flow, a match here implies that the data is being
        parsed accurately. A second species with no temperature
        intervals is also checked for, ensuring these special cases
        are also parsed as expected.
        
        """
        self.assertTrue(test_reactant in self.parsed['reactants'])
        self.assertTrue(test_sp_reactant in self.parsed['reactants'])


class TestLookup(unittest.TestCase):
    def test_lookup(self):
        """Test sample species is in lookup results"""
        self.assertTrue(test_reactant in 
                        thermoinp.lookup('JP-10')['reactants'])
    
    def test_lookup_exact(self):
        """Test a search for a single species with the 'exact' flag.
        
        This test checks that the exact flag correctly restricts
        results to only an exact name match, i.e. the 'gas_products'
        category is a list containing only the H2 dataset and not
        a list of datasets where 'H2' matches the start of the name
        (H2O2, etc). 

        """
        matches = thermoinp.lookup('H2', exact=True)['gas_products']
        self.assertEqual([test_gas], matches)

class TestCreateSubset(unittest.TestCase):
    def setUp(self):
        self.datad = os.path.join(os.path.dirname(__file__), 'data')

    def get_data(self, fname):
        """Return the contents of a specified ThermoBuild data file""" 
        path = os.path.join(self.datad, fname)
        with open(path, 'r') as f: return f.read().strip('\n')

    def test_no_reactants(self):
        """Test a subset containing only reactant/product species."""
        correct_data = self.get_data('products_subset.txt')
        species = ('O2', 'N2', 'Ar', 'CO2')
        testing_data = thermoinp.create_subset(species, exact=True)
        self.assertEqual(testing_data, correct_data)

    def test_no_products(self):
        """Test a subset containing only reactant species."""
        correct_data = self.get_data('reactants_subset.txt')
        species = ('Air', 'Jet-A(L)', 'Jet-A(g)')
        testing_data = thermoinp.create_subset(species, exact=True)
        self.assertEqual(testing_data, correct_data)
        
    def test_mixed(self):
        """Test a subset containing both reactants and products."""
        correct_data = self.get_data('mixed_subset.txt')
        species = ('C3H8', 'Air')
        testing_data = thermoinp.create_subset(species, exact=True)
        self.assertEqual(testing_data, correct_data)

# --------------------------------------------------------------------
# TEST DATA
# --------------------------------------------------------------------

test_gas = Species(
    name='H2',
    comments='Ref-Elm. Gurvich,1978 pt1 p103 pt2 p31.',
    nintervals=3,
    formula='H:2.00',
    phase=0,
    refcode='tpis78',
    molwt=2.0158800,
    h_formation=0.000,
    h_assigned=None,
    T_reference=None,
    intervals=[_Interval((200, 1000),
                         7,
                         (-2, -1, 0, 1, 2, 3, 4, 0),
                         8468.102,
                         ( 4.078323210e+04,
                          -8.009186040e+02,
                           8.214702010e+00,
                          -1.269714457e-02,
                           1.753605076e-05,
                          -1.202860270e-08,
                           3.368093490e-12,
                          ),
                         ( 2.682484665e+03,
                          -3.043788844e+01,
                          ),
                         ),
               _Interval((1000, 6000),
                         7,
                         (-2, -1, 0, 1, 2, 3, 4, 0),
                         8468.102,
                         ( 5.608128010e+05,
                          -8.371504740e+02,
                           2.975364532e+00,
                           1.252249124e-03,
                          -3.740716190e-07,
                           5.936625200e-11,
                          -3.606994100e-15,
                          ),
                         ( 5.339824410e+03,
                          -2.202774769e+00,
                          )
                         ),
               _Interval((6000.0,  20000.0),
                         7,
                         (-2, -1, 0, 1, 2, 3, 4, 0),
                         8468.102,
                         ( 4.966884120e+08,
                          -3.147547149e+05,
                           7.984121880e+01,
                          -8.414789210e-03,
                           4.753248350e-07,
                          -1.371873492e-11,
                           1.605461756e-16,
                           ),
                         ( 2.488433516e+06,
                          -6.695728110e+02,
                          )
                         )
                        ],
    )

test_condensed = Species(
    name='Ag(cr)',
    comments='Cubic. Ref-Elm. Cox,1989 p228.',
    nintervals=1,
    refcode='coda89',
    formula='AG:1.00',
    phase=1,
    molwt=107.8682000,
    h_formation=0.000,
    h_assigned=None,
    T_reference=None,
    intervals=[_Interval((200.0, 1235.080),
                         7,
                         (-2, -1, 0, 1, 2, 3, 4, 0),
                         5745.000,
                         (-7.099236470e+04,
                           7.254788020e+02,
                           1.066518380e-01,
                           5.529541550e-03,
                          -4.425590850e-06,
                           2.091668120e-09,
                          -3.888924460e-13,
                          ),
                         (-4.614014260e+03,
                           5.074216040e+00,
                          ),
                         ),  # end interval 1
                 ], # end intervals
    )

test_reactant = Species(
    name='JP-10(g)',
    comments=('Exo-tetrahydrodicyclopentadiene. Zehe,2002.'
             '             React.'),
    nintervals=2,
    refcode='g 6/01',
    formula='C:10.00 H:16.00',
    phase=0,
    molwt=136.2340400,
    h_formation=-86855.900,
    h_assigned=None,
    T_reference=None,
    intervals=[_Interval((200.0, 1000.0),
                         7,
                         (-2, -1, 0, 1, 2, 3, 4, 0),
                         22997.434,
                         (-7.310769440e+05,
                           1.521764245e+04,
                          -1.139312644e+02,
                           4.281501620e-01,
                          -5.218740440e-04,
                           3.357233400e-07,
                          -8.805750980e-11,
                          ),
                         (-8.067482120e+04,
                           6.320148610e+02,
                           )
                         ),  # end interval 1
               _Interval((1000.0, 6000.0),
                         7,
                         (-2, -1, 0, 1, 2, 3, 4, 0),
                         22997.434,
                         ( 1.220329594e+07,
                          -5.794846240e+04,
                           1.092281156e+02,
                          -1.082406215e-02,
                           2.034992622e-06,
                          -2.052060369e-10,
                           8.575760210e-15,
                           ),             
                         ( 3.257334050e+05,
                          -7.092350760e+02,
                           ),
                         ),  # end interval 2
               ],		  # end Intervals
    ) # end Species
                    

test_sp_reactant = Species(
    name='RP-1',
    comments=('Mehta et.al. AIAA 95-2962 1995.'
           ' Hcomb(high) = 19923.BTU/#'),
    nintervals=0,
    refcode='gll/00',
    formula='C:1.00 H:1.95',
    phase=1,
    molwt=13.9761830,
    h_formation=None,
    h_assigned=-24717.700,
    T_reference=298.150,
    intervals=None,
    )

if __name__ == '__main__':
    unittest.main()
