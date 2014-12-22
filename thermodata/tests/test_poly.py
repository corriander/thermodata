import unittest

from thermodata.poly import NASAPoly, NASAPolyND, NASAPolyML, Parser

# TODO: Fill out these tests for the NASAPoly variants.

class TestNASAPoly(unittest.TestCase):
    # There's nothing much to test here!

    def test_plot(self):
        self.skipTest("Functionality not implemented.")

    def test_linearise(self):
        """Could very well be possible to approximate into a linear
        equation for speeding up calculations."""
        self.skipTest("Functionality not implemented (may never be).")


class TestNASAPolyND(unittest.TestCase):

    def test_cpnd(self):
        # We want to do this for a range of different species/polys
        self.skipTest("Test not implemented.")


class TestNASAPolyML(unittest.TestCase):

    def test_cpmol(self):
        # We want to do this for a range of different species/polys
        self.skipTest("Test not implemented.")


class TestParser(unittest.TestCase):

    p = Parser()

    def test__double_array_to_float_single(self):
        strings = (
            "-1.008408053D-08",
            " 9.359095680D+05",
            "-5.525864490D+02",
        )
        floats = (
            -1.008408053e-8,
            9.359095680e5,
            -5.525864490e2,
        )

        for string, flt in zip(strings, floats):
            parsed = self.p._double_array_to_float(string)
            self.assertAlmostEqual(flt, *parsed)

    def test__double_array_to_float_record(self):

        record = gas2i.splitlines()[-1]
        floats = (
            -3.945148580e-11,
            1.509592679e-15,
            0.0,
            4.951216910e4,
            -5.417083590e1,
        )
        zobj = zip(self.p._double_array_to_float(record), floats)
        for a, b in zobj:
            self.assertAlmostEqual(a, b)

    def test__parse_coefficients(self):
        records = gas2i.splitlines()[-2:]

        coefficients = (
             9.359095680e05,
            -4.441073340e03,
             1.368958451e01,
            -1.647526929e-03,
             3.819863520e-07,
            -3.945148580e-11,
             1.509592679e-15,
        )
        constants = (
            4.951216910e4,
            -5.417083590e1,
        )
        a, b = self.p._parse_coefficients(records)

        for x, y in zip(a, coefficients):
            self.assertAlmostEqual(x, y)

        for x, y in zip(b, constants):
            self.assertAlmostEqual(x, y)

    def test__parse_metadata(self):
        record = gas2i.splitlines()[-3]

        lim, n, exp, dh = self.p._parse_metadata(record)
        self.assertEqual(lim, (1000, 6000))
        self.assertEqual(n, 7)
        self.assertEqual(exp, (-2, -1, 0, 1, 2, 3, 4, 0))
        self.assertAlmostEqual(dh, 13594.351)

    def test__parse_interval(self):
        records = gas2i.splitlines()[-3:]

        pobj = self.p._parse_interval(records)
        self.assertAlmostEqual(pobj.lim[1], 6000.)
        self.assertAlmostEqual(pobj.a[3], -1.647526929e-3)
        self.assertAlmostEqual(pobj.b[0], 4.951216910e4)
        self.assertEqual(pobj.n, 7)
        self.assertEqual(pobj.exp, (-2, -1, 0, 1, 2, 3, 4, 0))
        self.assertEqual(pobj.dh, 13594.351)

    def test__parse_interval_default_type(self):
        records = gas2i.splitlines()[-3:]

        pobj = self.p._parse_interval(records)
        self.assertIsInstance(pobj, NASAPoly)

    def test__parse_interval_nd_type(self):
        p = Parser(NASAPolyND)
        records = gas2i.splitlines()[-3:]

        pobj = p._parse_interval(records)
        self.assertIsInstance(pobj, NASAPolyND)

    # Don't bother testing intervals, it either works or it doesn't.

    def test__parse_intervals(self):
        records = gas2i.splitlines()[2:]
        self.assertEqual(len(self.p._parse_intervals(records)), 2)

    def test___call__gas2i(self):
        """Parser returns right num. polys for sample dataset."""
        pobjs = self.p(gas2i)
        self.assertEqual(len(pobjs), 2)

    # TODO: More tests for different datasets.


# --------------------------------------------------------------------
# Test Data
# --------------------------------------------------------------------
# 2-interval gas
gas2i = (
"""OCCN              Cyanooxomethyl radical. Dorofeeva,2001.
 2 srd 01 C   2.00N   1.00O   1.00    0.00    0.00 0   54.0275000     210000.000
    200.000  1000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13594.351
 2.428801128D+04-5.525864490D+02 8.587838090D+00-3.379387570D-03 1.119841795D-05
-1.008408053D-08 3.086448751D-12 0.000000000D+00 2.599620066D+04-1.659592428D+01
   1000.000  6000.000 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13594.351
 9.359095680D+05-4.441073340D+03 1.368958451D+01-1.647526929D-03 3.819863520D-07
-3.945148580D-11 1.509592679D-15 0.000000000D+00 4.951216910D+04-5.417083590D+01
""")

if __name__ == '__main__':
    unittest.main()
