import unittest

from thermodata.poly import NASAPoly, NASAPolyND, NASAPolyML

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


if __name__ == '__main__':
    unittest.main()
