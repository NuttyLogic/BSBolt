import subprocess
from typing import Dict, Union
from BSBolt.Utils.UtilityFunctions import retrieve_iupac


class StreamSim:

    def __init__(self, paired_end=False, sim_command=None):
        self.paired_end = paired_end
        self.sim_command = sim_command
        self.contig_variants = {}
        self.variant_contig = None

    def __iter__(self):
        sim = subprocess.Popen(self.sim_command,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
        variant_output = False
        sim_output = iter(sim.stdout.readline, '')
        read_pair = {1: None, 2: None}
        paired_count = 0
        while True:
            try:
                formatted_line = next(sim_output).strip()
            except StopIteration:
                break
            else:
                if formatted_line == 'Contig Variant Start':
                    variant_output = True
                elif variant_output:
                    variant_output = self.collect_variant_info(formatted_line)
                    if not variant_output and self.contig_variants:
                        yield self.variant_contig, self.contig_variants
                        self.contig_variants = {}
                else:
                    read_info = self.process_read_name(formatted_line)
                    seq = next(sim_output).strip()
                    comment = next(sim_output).strip()
                    qual = self.modify_qual(next(sim_output).strip())
                    read_info.update(dict(comment=comment, seq=seq, qual=qual))
                    read_pair[read_info['pair']] = read_info
                    paired_count += 1
                    if paired_count == 2:
                        assert read_pair[1]['read_id'] == read_pair[2]['read_id']
                        yield False, read_pair
                        read_pair = {1: None, 2: None}
                        paired_count = 0

    @staticmethod
    def modify_qual(quality: str) -> str:
        # modify start position to ensure proper qual handling downstream
        quality_split = list(quality)
        qual_start = int(ord(quality_split[0]) - 33)
        quality_split[0] = chr(qual_start + 32)
        return ''.join(quality_split)

    def collect_variant_info(self, formatted_line):
        if formatted_line == 'Contig Variant End':
            return False
        variant_info = self.process_variant_line(formatted_line)
        assert variant_info['pos'] not in self.contig_variants
        self.contig_variants[variant_info['pos']] = variant_info
        self.variant_contig = variant_info['chrom']
        return True

    @staticmethod
    def process_variant_line(formatted_line):
        chrom, pos, reference, alt, heterozygous = formatted_line.split('\t')
        pos = int(pos)
        indel = 0
        iupac = None
        if reference == '-':
            indel = 1
        elif alt == '-':
            indel = -1
        else:
            iupac = retrieve_iupac(alt)
        return dict(chrom=chrom, pos=pos, reference=reference, alt=alt, heterozygous=heterozygous,
                    indel=indel, iupac=iupac)

    @staticmethod
    def process_read_name(formatted_line: str) -> Dict[str, Union[str, int]]:
        read_info = formatted_line.split(':')
        chrom, start, end, insert_size, read_id, cigar, pair, c_base_info, g_base_info = read_info
        return dict(chrom=chrom.replace('@', ''), start=int(start), end=int(end),
                    insert_size=insert_size, read_id=read_id, cigar=cigar, pair=int(pair),
                    c_base_info=c_base_info, g_base_info=g_base_info)
