from collections import Counter
import pickle
import pysam


class CallMethylationValues:
    """
    Keyword Arguments
        input_file (str): str to input bam/sam file
        genome_database (str): str to genome directory
        rm_overlap (bool):  remove overlapping signal, by default calls from one overlapping read
        remove_ccgg (bool): don't call CCGG sequences
        max_read_depth (int): default = 8000
        contig (str): contig to process
        min_base_quality (int): minimum quality for a base to be reported for methylation calling
    Attributes:
        self.input_file (str): path to input bam/sam file
        self.input_bam (pysam.Samfile): pysam.Samfile object to retrieve pileup information
        self.genome_database (str): formatted path to genome database
        self.rm_overlap (bool):  remove overlapping signal, by default calls from one overlapping read
        self.remove_ccgg (bool): don't call CCGG sequences
        self.max_read_depth (int): default = 8000
        self.contig (str): contig to process
        self.min_base_quality (int): minimum quality for a base to be reported for methylation calling
        self.context_tables (dict): dict of dicts listing nucleotide context
    """

    def __init__(self, input_file: str = None, genome_database: str = None,
                 ignore_overlap: bool = True, remove_ccgg: bool = False, ignore_orphans: bool = True,
                 max_read_depth: int = 8000, contig: str = None, min_base_quality: int = 1, return_queue=None,
                 cg_only: bool = False):
        self.input_file = str(input_file)
        self.input_bam = pysam.AlignmentFile(self.input_file, mode='rb',
                                             require_index=True)
        self.genome_database = str(genome_database)
        if self.genome_database[-1] != '/':
            self.genome_database = f'{self.genome_database}/'
        self.ignore_overlap = ignore_overlap
        self.remove_ccgg = remove_ccgg
        self.ignore_orphans = ignore_orphans
        self.max_read_depth = max_read_depth
        self.contig = contig
        self.min_base_quality = min_base_quality
        self.chunk_size = 10000
        self.cg_only = cg_only
        self.context_tables = self.get_context_tables
        self.return_queue = return_queue
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
        try:
            chrom_seq = self.get_reference_sequence(f'{self.genome_database}{self.contig}.pkl')
        except FileNotFoundError:
            print(f'{self.contig} not found in BSBolt DB, Methylation Calls for {self.contig} skipped. Methylation '
                  f'values should be called using the same DB used for alignment.')
            self.return_queue.put([])
        else:
            self.call_contig(chrom_seq)

    def call_contig(self, chrom_seq):
        """Iterates through bam pileup, calling methylation values if the reference nucleotide is a C or G. Pileup reads
        are buffered and accessed as needed.
        """
        # iterate through pileup
        line_count = 0
        contig_chunk = []
        # flag_require (0) mapped read
        # flag_filter 1540 = pcr duplicates (1024) + read unmapped (4) + read fails platform/vendor quality checks (512)
        for pileup_col in self.input_bam.pileup(max_depth=self.max_read_depth,
                                                contig=self.contig,
                                                ignore_overlaps=self.ignore_overlap,
                                                min_base_quality=self.min_base_quality,
                                                ignore_orphans=self.ignore_orphans,
                                                flag_require=0,
                                                flag_filter=1540):
            # get sequence around pileup site
            reference_seq = chrom_seq[(pileup_col.reference_pos - 3):(pileup_col.reference_pos + 4)].upper()
            # get nucleotide context
            fivemer = reference_seq[1:-1]
            if len(fivemer) == 5:
                nucleotide = fivemer[2]
                context, subcontext = self.get_context(nucleotide, fivemer)
                if self.cg_only and subcontext != 'CG':
                    continue
                # check if CCGG in sequence, skip loop if filter True
                if self.check_ccgg(reference_seq):
                    continue
                # count the pileup read bases, Uppercase watson strand, lowercase crick strand
                try:
                    base_counts = Counter(pileup_col.get_query_sequences(mark_matches=False,
                                                                         mark_ends=False,
                                                                         add_indels=False))
                except AssertionError:
                    # pysam may throw an error if the number of reads in pileup column is above the max
                    continue
                meth_line = self.get_methylation_call(nucleotide, base_counts)
                meth_line.update({'pos': pileup_col.reference_pos + 1, 'chrom': self.contig,
                                  'context': context, 'subcontext': subcontext})
                contig_chunk.append(tuple(meth_line.values()))
                line_count += 1
                # return chunk for output
                if line_count == self.chunk_size - 1:
                    self.return_queue.put(contig_chunk, block=True)
                    contig_chunk = []
                    line_count = 0
        self.return_queue.put(contig_chunk, block=True)

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

    def check_ccgg(self, sequence):
        """checks if sequence is == to CCGG
        """
        if self.remove_ccgg:
            return 'CCGG' in sequence
        return False

    @staticmethod
    def get_methylation_call(nucleotide, base_counts):
        """
        Methylation for a C relative to the sense strand of the reference can only be called using watson reads,
        and G with crick reads
        Arguments
            nucleotide (str): reference nucleotide
            base_counts (collections.Counter): watson nucleotides are Uppercase and crick nucleotides lowercase
        Returns:
             methylation call dictionary
        """
        meth_cytosines = 0
        unmeth_cytosines = 0
        # call cytonsines with watson
        if nucleotide == 'C':
            meth_cytosines = base_counts.get('C', 0)
            unmeth_cytosines = base_counts.get('T', 0)
        # call guanines with crick strand
        elif nucleotide == 'G':
            meth_cytosines = base_counts.get('g', 0)
            unmeth_cytosines = base_counts.get('a', 0)
        all_cytosines = meth_cytosines + unmeth_cytosines
        meth_level = 'na'
        if all_cytosines > 0:
            # if cytosines present call methylation
            meth_level = round(float(meth_cytosines) / float(all_cytosines), 3)
        f'{base_counts.get("A", 0)}\t{base_counts.get("T", 0)}\t{base_counts.get("C", 0)}' \
            f'\t{base_counts.get("G", 0)}\t{base_counts.get("N", 0)}'
        forward_counts = f'{base_counts.get("A", 0)}\t{base_counts.get("T", 0)}\t{base_counts.get("C", 0)}' \
            f'\t{base_counts.get("G", 0)}\t{base_counts.get("N", 0)}'
        reverse_counts = f'{base_counts.get("a", 0)}\t{base_counts.get("t", 0)}\t{base_counts.get("c", 0)}' \
            f'\t{base_counts.get("g", 0)}\t{base_counts.get("n", 0)}'
        return {'nucleotide': nucleotide, 'meth_cytosines': meth_cytosines, 'unmeth_cytosines': unmeth_cytosines,
                'all_cytosines': all_cytosines, 'meth_level': meth_level, 'forward_counts': forward_counts,
                'reverse_counts': reverse_counts}

    @staticmethod
    def get_reference_sequence(path):
        """load serialized reference file from path
        """
        with open(path, 'rb') as genome_file:
            return pickle.load(genome_file)
