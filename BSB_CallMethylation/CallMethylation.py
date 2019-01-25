import pickle
import pysam
import time


def time_function(start, message):
    total_time = time.time() - start
    print(f'{message} {total_time}')
    return time.time()


class CallMethylation:
    """
    Keyword Arguments
        input_file (str): str to input bam/sam file
        genome_database (str): str to genome directory
        remove_sx_reads (bool): remove incompletely converted reads
        rm_overlap (bool):  remove overlapping signal, by default calls from one overlapping read
        remove_ccgg (bool): don't call CCGG sequences
        min_read_depth (int): default = 1
        max_read_depth (int): default = 8000
        contig (str): contig to process
        min_base_quality (int): minimum quality for a base to be reported for methylation calling
    Attributes:
        self.input_file (str): path to input bam/sam file
        self.input_bam (pysam.Samfile): pysam.Samfile object to retrieve pileup information
        self.genome_database (str): formatted path to genome database
        self.remove_sx_reads (bool): remove incompletely converted reads
        self.rm_overlap (bool):  remove overlapping signal, by default calls from one overlapping read
        self.remove_ccgg (bool): don't call CCGG sequences
        self.min_read_depth (int): default = 1
        self.max_read_depth (int): default = 8000
        self.contig (str): contig to process
        self.min_base_quality (int): minimum quality for a base to be reported for methylation calling
        self.context_tables (dict): dict of dicts listing nucleotide context
    """

    def __init__(self, input_file=None, genome_database=None,
                 remove_sx_reads=True, ignore_overlap=False, remove_ccgg=False,
                 min_read_depth=10, max_read_depth=8000, contig=None, min_base_quality=0):
        assert isinstance(input_file, str), 'Path to input file not valid'
        assert isinstance(genome_database, str), 'Path to genome database not valid'
        assert isinstance(remove_sx_reads, bool), 'Not valid bool'
        assert isinstance(ignore_overlap, bool), 'Not valid bool'
        assert isinstance(remove_ccgg, bool), 'Not valid bool'
        assert isinstance(min_read_depth, int), 'Minimum depth must be integer'
        assert isinstance(max_read_depth, int), 'Maximum depth must be integer'
        assert isinstance(contig, str), 'Contig must be string'
        self.input_file = str(input_file)
        self.input_bam = pysam.AlignmentFile(self.input_file, 'rb')
        self.genome_database = str(genome_database)
        if self.genome_database[-1] != '/':
            self.genome_database = f'{self.genome_database}/'
        self.remove_sx_reads = remove_sx_reads
        self.ignore_overlap = ignore_overlap
        self.remove_ccgg = remove_ccgg
        self.min_read_depth = min_read_depth
        self.max_read_depth = max_read_depth
        self.contig = contig
        self.min_base_quality = min_base_quality
        self.context_tables = self.get_context_tables
        self.return_list = []
        self.counting_dict = {}

    @property
    def get_context_tables(self):
        """
        Returns:
            context_tables
        """
        context_tables = {'context_table': {'CAA': 'CHH', 'CAC': 'CHH', 'CAG': 'CHG', 'CAT': 'CHH', 'CCA': 'CHH',
                                            'CCC': 'CHH', 'CCG': 'CHG', 'CCT': 'CHH', 'CGA': 'CG', 'CGC': 'CG',
                                            'CGG': 'CG', 'CGT': 'CG', 'CTA': 'CHH', 'CTC': 'CHH', 'CTG': 'CHG',
                                            'CTT': 'CHH'},
                          'sub_context_table': {'CAA': 'CA', 'CAC': 'CA', 'CAG': 'CA', 'CAT': 'CA', 'CCA': 'CC',
                                                'CCC': 'CC', 'CCG': 'CC', 'CCT': 'CC', 'CGA': 'CG', 'CGC': 'CG',
                                                'CGG': 'CG', 'CGT': 'CG', 'CTA': 'CT', 'CTC': 'CT', 'CTG': 'CT',
                                                'CTT': 'CT'},
                          'antisense_context_table': {'TTG': 'CHH', 'TGG': 'CHH', 'TCG': 'CG', 'TAG': 'CHH',
                                                      'GTG': 'CHH', 'GGG': 'CHH', 'GCG': 'CG', 'GAG': 'CHH',
                                                      'CTG': 'CHG', 'CGG': 'CHG', 'CCG': 'CG', 'CAG': 'CHG',
                                                      'ATG': 'CHH', 'AGG': 'CHH', 'ACG': 'CG', 'AAG': 'CHH'},
                          'antisense_sub_context_table': {'TTG': 'CA', 'TGG': 'CC', 'TCG': 'CG', 'TAG': 'CT',
                                                          'GTG': 'CA', 'GGG': 'CC', 'GCG': 'CG', 'GAG': 'CT',
                                                          'CTG': 'CA', 'CGG': 'CC', 'CCG': 'CG', 'CAG': 'CT',
                                                          'ATG': 'CA', 'AGG': 'CC', 'ACG': 'CG', 'AAG': 'CT'}}
        return context_tables

    def call_methylation(self):
        """Iterates through bam pileup, calling methylation values if the reference nucleotide is a C or G. Pileup reads
        are bufferend and accessed as needed then deleted when they exit the scope.
        @return:
        """
        # load serialized reference sequence
        chrom_seq = self.get_reference_sequence(f'{self.genome_database}{self.contig}.pkl')
        # iterate through pileup
        for pileup_col in self.input_bam.pileup(max_depth=self.max_read_depth, 
                                                contig=self.contig,
                                                ignore_overlaps=self.ignore_overlap, 
                                                min_base_quality=self.min_base_quality):
            # initialize count dictionaries for pileup sites
            ATCG_forward = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
            ATCG_reverse = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
            # get sequence around pileup site
            reference_seq = chrom_seq[(pileup_col.reference_pos - 3):(pileup_col.reference_pos + 4)].upper()
            # get nucleotide context
            fivemer = reference_seq[1:-1]
            if len(fivemer) == 5:
                nucleotide = fivemer[2]
                context, subcontext = self.get_context(nucleotide, fivemer)

                # check if sequence is CCGG, skip loop if filter True
                if self.check_ccgg(reference_seq):
                    continue

                for pileup_read in pileup_col.pileups:
                    if not self.check_read(pileup_read):
                        continue
                    if pileup_read.query_position:
                        read_nucleotide = pileup_read.alignment.query_sequence[pileup_read.query_position]
                        if pileup_read.alignment.is_reverse:
                            ATCG_reverse[read_nucleotide] += 1
                        else:
                            ATCG_forward[read_nucleotide] += 1
                # get methylation call at line
                meth_line = self.get_methylation_call(nucleotide, ATCG_forward, ATCG_reverse)
                meth_line.update({'pos': pileup_col.reference_pos + 1, 'chrom': self.contig,
                                  'context': context, 'subcontext': subcontext})
                self.return_list.append(meth_line)

    def check_read(self, pileup_read):
        """Check if read converted and pileup_read location and indel
        """
        # proceed after conversion check
        if self.check_conversion(pileup_read) and not pileup_read.indel:
            return True
        return False

    def get_context(self, nucleotide, fivemer):
        """
        Arguments:
            nucleotide (str): 1 nucleotide
            fivemer (str): 5 nucleotides
        Returns:
            context (str): methylation context
            subcontext (str): methylation subcontext
        """
        null_context = '--'
        context = null_context
        subcontext = null_context
        if nucleotide == 'C':
            subcontext = self.context_tables['sub_context_table'].get(fivemer[2:5], null_context)
            context = self.context_tables['context_table'].get(fivemer[2:5], null_context)
        elif nucleotide == 'G':
            subcontext = self.context_tables['antisense_sub_context_table'].get(fivemer[0:3], null_context)
            context = self.context_tables['antisense_context_table'].get(fivemer[0:3], null_context)
        return context, subcontext

    def check_conversion(self, pileup_read):
        """Checks if XS flag present
        """
        if self.remove_sx_reads:
            return ('XS', 1) not in pileup_read.alignment.tags
        return True

    def check_ccgg(self, sequence):
        """checks if sequence is == to CCGG
        """
        if self.remove_ccgg:
            return 'CCGG' in sequence
        return False

    @staticmethod
    def get_methylation_call(nucleotide, forward_meth_dict, reverse_meth_dict):
        """
        Arguments
            nucleotide (str): reference nucleotide
            forward_meth_dict (dict): dictionary of bases called from sense strand
            reverse_meth_dict (dict): dictionary of bases called from antisense strand
        Returns:
             methylation call dictionary
        """
        meth_cytosines = 0
        unmeth_cytosines = 0
        if nucleotide == 'C':
            meth_cytosines = forward_meth_dict['C']
            unmeth_cytosines = forward_meth_dict['T']
        elif nucleotide == 'G':
            meth_cytosines = reverse_meth_dict['G']
            unmeth_cytosines = reverse_meth_dict['A']
        all_cytosines = meth_cytosines + unmeth_cytosines
        meth_level = 'na'
        if all_cytosines > 0:
            # if cytosines present call methylation
            meth_level = round(float(meth_cytosines) / float(all_cytosines), 3)
        forward_counts = '\t'.join([str(count) for count in forward_meth_dict.values()])
        reverse_counts = '\t'.join([str(count) for count in reverse_meth_dict.values()])
        return {'nucleotide': nucleotide, 'meth_cytosines': meth_cytosines, 'unmeth_cytosines': unmeth_cytosines,
                'all_cytosines': all_cytosines, 'meth_level': meth_level, 'forward_counts': forward_counts,
                'reverse_counts': reverse_counts}

    @staticmethod
    def get_reference_sequence(path):
        """load serialized reference file from path
        """
        with open(path, 'rb') as genome_file:
            return pickle.load(genome_file)
