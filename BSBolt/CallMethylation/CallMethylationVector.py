import pickle
import pysam
import numpy as np


class CallMethylationVector:
    """
     Experimental
     return vector of methylated sites
    """

    def __init__(self, input_file: str = None, genome_database: str = None,
                 max_read_depth: int = 10000, contig: str = None, min_base_quality: int = 0, return_queue=None,
                 cg_only: bool = False, start=None, end=None):
        self.input_file = str(input_file)
        self.input_bam = pysam.AlignmentFile(self.input_file, 'rb')
        self.genome_database = str(genome_database)
        if self.genome_database[-1] != '/':
            self.genome_database = f'{self.genome_database}/'
        self.max_read_depth = max_read_depth
        self.contig = contig
        self.start = start
        self.end = end
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
            self.return_queue.put([])
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
        methylation_vectors = {}
        for pileup_col in self.input_bam.pileup(max_depth=self.max_read_depth,
                                                contig=self.contig,
                                                start=self.start,
                                                end=self.end,
                                                ignore_overlaps=False,
                                                min_base_quality=self.min_base_quality):
            # get sequence around pileup site
            reference_seq = chrom_seq[(pileup_col.reference_pos - 3):(pileup_col.reference_pos + 4)].upper()
            # get nucleotide context
            fivemer = reference_seq[1:-1]
            if len(fivemer) == 5:
                nucleotide = fivemer[2]
                context, subcontext = self.get_context(nucleotide, fivemer)

                if self.cg_only and subcontext != 'CG':
                    continue
                # count the pileup read bases, Uppercase watson strand, lowercase crick strand
                try:
                    base_counts = pileup_col.get_query_sequences(mark_matches=False,
                                                                 mark_ends=False,
                                                                 add_indels=False)
                    read_names = pileup_col.get_query_names()
                except AssertionError:
                    # pysam may throw an error if the number of reads in pileup column is above the max
                    continue
                meth_calls = self.get_methylation_call(nucleotide, read_names, base_counts)
                for read, call in meth_calls:
                    if read in methylation_vectors:
                        methylation_vectors[read]['pos'].append(pileup_col.reference_pos)
                        methylation_vectors[read]['calls'].append(call)
                    else:
                        methylation_vectors[read] = {'pos': [pileup_col.reference_pos],
                                                     'calls': [call],
                                                     'chrom': self.contig}
        for read in methylation_vectors.keys():
            methylation_vectors[read]['pos'] = np.array(methylation_vectors[read]['pos'])
            methylation_vectors[read]['calls'] = np.array(methylation_vectors[read]['calls'])
            methylation_vectors[read]['range'] = np.array([methylation_vectors[read]['pos'][0],
                                                           methylation_vectors[read]['pos'][-1]])
        self.return_queue.append(methylation_vectors)

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

    @staticmethod
    def get_methylation_call(nucleotide, read_name, base_calls):
        """
        Methylation for a C relative to the sense strand of the reference can only be called using watson reads,
        and G with crick reads
        Arguments
            nucleotide (str): reference nucleotide
            base_counts (collections.Counter): watson nucleotides are Uppercase and crick nucleotides lowercase
        Returns:
             methylation call dictionary
        """
        # call cytonsines with watson
        called_reads = []
        for read, call in zip(read_name, base_calls):
            if nucleotide == 'C':
                if call == 'C' or call == 'G':
                    called_reads.append((read, 1))
                elif call == 'T' or call == 'A':
                    called_reads.append((read, 0))
            elif nucleotide == 'G':
                if call == 'c' or call == 'g':
                    called_reads.append((read, 1))
                elif call == 't' or call == 'a':
                    called_reads.append((read, 0))
        return called_reads

    @staticmethod
    def get_reference_sequence(path):
        """load serialized reference file from path
        """
        with open(path, 'rb') as genome_file:
            return pickle.load(genome_file)
