from BSBolt.Utils.UtilityFunctions import reverse_complement, retrieve_iupac


class ProcessCutSites:
    """ Process cut format information and returns a dictionary with restriction sequence, reverse complement of
    restriction sequence and offset.
    Keyword Arguments:
        cut_format (str):
        cut_descriptor (str):
    Attributes:
        self.cut_format (list): list of recognition sites split by ,
        self.cut_descriptor (str): default = -, character that specifies where the restriction enzyme cuts in the
            recognition sequence
        self.restriction_site_dict (dict): dictionary of restriction sites with offset information
        self.process_cut_sites (func): population self.restriction_site_dict

    """

    def __init__(self, cut_format=None, cut_descriptor='-'):
        assert isinstance(cut_format, str), 'cut_format must be a str'
        self.cut_format: list = cut_format.upper().replace(' ', '').split(',')
        self.cut_descriptor = cut_descriptor
        self.restriction_site_dict = {}
        self.process_cut_sites()

    def process_cut_sites(self):
        """Process forward and reverse complement of all restriction sites
        """
        for site in self.cut_format:
            # get offset for restriction sequence
            forward_offset, reverse_offset = self.get_site_offsets(site)
            # remove cut descriptor from recognition sequence
            nucleotide_site = site.replace(self.cut_descriptor, '')
            # if iupac symbol in restriction site get all possible recognition sequences
            forward_sites, reverse_sites = self.get_recognition_site_sequences(nucleotide_site)
            # iterate through site and population self.restriction_site_dict
            for f_site, r_site in zip(forward_sites, reverse_sites):
                self.restriction_site_dict[f_site] = forward_offset
                # if the f_site is equivalent to the r_site don't add reverse_site as the offset will be incorrect
                if f_site != r_site:
                    self.restriction_site_dict[r_site] = reverse_offset

    def get_site_offsets(self, site):
        """Given a restriction site return offset, or proper alignment based on where cut occurs in sequence
        Arguments:
            site (str): recognition site sequence
        Returns:
            forward_offest (int/bool): if cut descriptor int, else False
            reverse_offset (int/bool): if cut descriptor int, else False
        """
        try:
            # retrieve index of cut descriptor where present
            forward_offset = site.index(self.cut_descriptor)
        except ValueError:
            # if descriptor not present return False
            forward_offset = False
            reverse_offset = False
        else:
            # get location of cut_descriptor in reverse complement
            reverse_offset = reverse_complement(site).index(self.cut_descriptor)
        return forward_offset, reverse_offset

    def get_recognition_site_sequences(self, recognition_site):
        """
        Arguments:
            recognition_site (str): sequence of site
        Returns:
            forward_recognition_sites (list): list of str
            reverse_recognition_sites (list): list of str
        """
        forward_recognition_sites = ['']
        # iterate through nucleotides in recognition site and replace cut_descriptor
        for nucelotide in recognition_site.replace(self.cut_descriptor, ''):
            temp_sites = []
            # retrieve_iupac returns value of iupac identifiers or the original nucleotide
            for iupac_nucleotide in retrieve_iupac(nucelotide):
                for site in forward_recognition_sites:
                    temp_sites.append(f'{site}{iupac_nucleotide}')
            forward_recognition_sites = temp_sites
        # get reverse complement for all forward recognition sites
        reverse_recognition_sites = [reverse_complement(site) for site in forward_recognition_sites]
        return forward_recognition_sites, reverse_recognition_sites
