import subprocess
from typing import Dict, Union
from BSBolt.Utils.UtilityFunctions import retrieve_iupac


class StreamSim:

    def __init__(self, paired_end=False):
        self.paired_end = paired_end
        self.contig_variants = {}

    def stream_sim_reads(self, sim_command):
        sim = subprocess.Popen(sim_command,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
        variant_output = False
        sim_output = iter(sim.stdout.readline, '')
        print(sim_command)
        print('here')
        while True:
            formatted_line = next(sim_output).strip()
            if not formatted_line:
                break
            if formatted_line == 'Contig Variant Start':
                variant_output = True
            elif variant_output:
                variant_output = self.collect_variant_info(formatted_line)
                if not variant_output:
                    yield False, self.contig_variants
                    self.contig_variants = {}
            else:
                read_info = self.process_read_name(formatted_line)
                seq = next(sim_output).strip()
                comment = next(sim_output).strip()
                qual = next(sim_output).strip()
                read_info.update(dict(comment=comment, seq=seq, qual=qual))
                if not self.paired_end and read_info['pair'] == 2:
                    continue
                yield True, read_info

    def collect_variant_info(self, formatted_line):
        if formatted_line == 'Contig Variant End':
            return False
        variant_info = self.process_variant_line(formatted_line)
        variant_label = f'{variant_info[0]}:{variant_info[1]}'
        assert variant_label not in self.contig_variants
        self.contig_variants[variant_label] = variant_info
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
        return chrom, pos, reference, alt, heterozygous, indel, iupac

    @staticmethod
    def process_read_name(formatted_line: str) -> Dict[str, Union[str, int]]:
        read_info = formatted_line.split(':')
        chrom, start, end, err_1, subs_1, indel_1, err_2, subs_2, indel_2, strand, read_id, pair = read_info
        return dict(chrom=chrom, start=int(start), end=int(end),
                    err_1=int(err_1), subs_1=int(subs_1), indel_1=int(indel_1),
                    err_2=int(err_2), subs_2=int(subs_2), indel_2=int(indel_2),
                    strand=int(strand), read_id=read_id, pair=int(pair))


'''
1)  contig name (chromsome name)
                2)  start end 1 (one-based)
                3)  end end 2 (one-based)
                4)  number of errors end 1
                5)  number of substitutions end 1
                6)  number of indels end 1
                5)  number of errors end 2
                6)  number of substitutions end 2
                7)  number of indels end 2
                8) strand 
                9) id
                10) pair
                '''