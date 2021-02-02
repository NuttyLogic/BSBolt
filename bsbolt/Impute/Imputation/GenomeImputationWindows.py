#! /user/bin/env python3

from typing import List, Tuple


def chrom_site_split(chrom_site: str) -> Tuple[str, int]:
    site_split = chrom_site.split(':')
    assert len(site_split) == 2, 'Site must be formatted as chromosome:pos'
    try:
        position = int(site_split[1])
    except ValueError as e:
        print('Genome position must be numeric')
        raise e
    return site_split[0], position


class GenomeImputationWindows:
    """Calculates overlapping windows for imputation, assumes sites are sorted by coordinate
        Keyword Arguments:
            site_labels (list): list of genomic sites, chr:pos
            imputation_window_size (int): size of pairwise distance window
        Attributes:
            self.imputation_ws (int): size of window to impute values (middle third of window)
            self.imputation_boundary (int): size of window step size when getting overlapping windows
            self.site_dictionary (dict): sliding window stats
            self.window_count (int): increasing count of imputation windows, doesn't include global window
            self.windows (list): list of site indices for windows
            self.site_window_dict (dict): site hashed to imputation window
            self.rolling_windows (dict): site for rolling windows
        """

    def __init__(self, site_labels: List[str] = None, imputation_window_size: int = 3000000):
        self.site_labels = site_labels
        self.imputation_ws = int(imputation_window_size / 3) * 3
        self.imputation_boundary = int(imputation_window_size / 3)
        self.site_dictionary = None
        self.boundary_trigger = False
        self.window_count = 0
        self.windows = []
        self.site_window_dict = {}
        self.rolling_windows = {'window1': [], 'window2': [], 'window3': []}
        self.run()

    def run(self):
        self.set_initial_state()
        self.get_windows()

    def set_initial_state(self):
        """Set imputation starting point"""
        chrom, site = chrom_site_split(self.site_labels[0])
        self.site_dictionary = {'Chromosome': chrom, 'SitePosition': site,
                                'Window1Boundary': site + self.imputation_boundary,
                                'Window2Boundary': site + 2 * self.imputation_boundary,
                                'Window3Boundary': site + self.imputation_ws}

    def site_window(self, site):
        position = '%s:%s' % (self.site_dictionary['Chromosome'], str(site))
        if site <= self.site_dictionary['Window1Boundary']:
            self.rolling_windows['window1'].append(position)
            self.site_window_dict[position] = self.window_count
        elif self.site_dictionary['Window1Boundary'] < site <= self.site_dictionary['Window2Boundary']:
            self.rolling_windows['window1'].append(position)
            self.rolling_windows['window2'].append(position)
            self.site_window_dict[position] = self.window_count
        elif self.site_dictionary['Window2Boundary'] < site <= self.site_dictionary['Window3Boundary']:
            self.rolling_windows['window1'].append(position)
            self.rolling_windows['window2'].append(position)
            self.rolling_windows['window3'].append(position)
            self.site_window_dict[position] = self.window_count + 1
        elif site > self.site_dictionary['Window3Boundary']:
            self.boundary_trigger = True

    def advance_windows(self):
        self.window_count += 1
        self.windows.append(self.rolling_windows['window1'])
        self.rolling_windows['window1'] = list(self.rolling_windows['window2'])
        self.rolling_windows['window2'] = list(self.rolling_windows['window3'])
        self.rolling_windows['window3'] = []
        self.site_dictionary['Window1Boundary'] = int(self.site_dictionary['Window2Boundary'])
        self.site_dictionary['Window2Boundary'] = int(self.site_dictionary['Window3Boundary'])
        self.site_dictionary['Window3Boundary'] = int(self.site_dictionary['Window3Boundary']
                                                      + self.imputation_boundary)

    def reset_windows(self, chrom, site):
        """Save current window and reset window at new position"""
        self.windows.append(self.rolling_windows['window1'])
        self.window_count += 1
        self.site_dictionary = {'Chromosome': chrom, 'SitePosition': site,
                                'Window1Boundary': site + self.imputation_boundary,
                                'Window2Boundary': site + 2 * self.imputation_boundary,
                                'Window3Boundary': site + self.imputation_ws}
        self.rolling_windows = {'window1': [], 'window2': [], 'window3': []}
        position = '%s:%s' % (chrom, str(site))
        self.rolling_windows['window1'].append(position)
        self.site_window_dict[position] = self.window_count

    def get_windows(self):
        """Rules for window extension and window reset"""
        for position in self.site_labels:
            chrom, site = chrom_site_split(position)
            if chrom != self.site_dictionary['Chromosome']:
                self.reset_windows(chrom, site)
            else:
                # raise assertion error if file not sorted
                assert site >= self.site_dictionary['SitePosition'], f'File must be position sorted; ' \
                                                                     f'{chrom}:{site} < ' \
                                                                     f'{chrom}:{self.site_dictionary["SitePosition"]}'
                self.site_window(site)
                if self.boundary_trigger:
                    # advance windows
                    self.advance_windows()
                    self.boundary_trigger = False
                    self.site_window(site)
                    if self.boundary_trigger:
                        self.reset_windows(chrom, site)
                        self.boundary_trigger = False
        self.windows.append(self.rolling_windows['window1'])
        self.windows.append(self.rolling_windows['window2'])
