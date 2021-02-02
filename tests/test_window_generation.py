from bsbolt.Impute.Imputation.GenomeImputationWindows import GenomeImputationWindows, chrom_site_split
import unittest

test_labels = [f'chr1:{x + 1}' for x in range(1000)]
test_labels.append('chr1:2000')
test_labels.extend(f'chr2:{x + 1}' for x in range(1000))

test_window = GenomeImputationWindows(site_labels=test_labels, imputation_window_size=300)


class TestSlidingWindow(unittest.TestCase):

    def setUp(self):
        pass

    def test_window_membership(self):
        self.assertIn(test_window.windows[0][-1], test_window.windows[1])
        self.assertIn(test_window.windows[0][-1], test_window.windows[2])
        self.assertNotIn(test_window.windows[0][-1], test_window.windows[3])

    def test_window_size(self):
        for window in test_window.windows:
            window_start = int(window[0].split(':')[1])
            window_end = int(window[-1].split(':')[1])
            self.assertLessEqual(window_end - window_start, 300)

    def test_window_locations(self):
        self.assertLessEqual(len(test_window.windows[9]), 1)

    def test_sort_order(self):
        try:
            GenomeImputationWindows(site_labels=['chr1:10', 'chr1:9'], imputation_window_size=300)
        except AssertionError as e:
            self.assertEqual(str(e), 'File must be position sorted; chr1:9 < chr1:10')

    def test_chromosome_membership(self):
        for window in test_window.windows:
            chrom = window[0].split(':')[0]
            for site in window:
                self.assertIn(chrom, site)

    def test_site_splitting(self):
        try:
            chrom_site_split('chr1000')
        except AssertionError as e:
            self.assertEqual(str(e), 'Site must be formatted as chromosome:pos')

    def test_numeric_position(self):
        test_value = 0
        try:
            chrom_site_split('chr:abc')
        except ValueError as e:
            self.assertEqual(str(e), 'invalid literal for int() with base 10: \'abc\'')
            test_value = 1
        self.assertEqual(test_value, 1)


if __name__ == '__main__':
    unittest.main()
