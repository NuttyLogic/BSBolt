#! /usr/env python3


class OpenSam:
    """Class to iterate through sam files and yield pairs of lines for paired end data and individual line for
    single end data
    Keyword Arguments:
        sam_file (str): path to sam file
        paired_end (bool):
    Attributes:
        self.f (TextIO): TextIO instance of sam file
        self.paired_end (bool): paired end data
    """

    def __init__(self, sam_file=None, paired_end=False):
        self.f = open(sam_file, 'r')
        self.paired_end = paired_end

    def __iter__(self):
        with self.f as sam:
            # keep current sam_id, if the same sam_id is encountered again don't yield this suppresses secondary
            # alignments, which is necessary to keep iteration order constant, possibly update in future versions
            sam_id = None
            while True:
                sam_line1 = sam.readline()
                # break if empty line encountered
                if not sam_line1:
                    break
                sam_dict1 = self.process_sam_line(sam_line1)
                if self.paired_end:
                    # reads should always be paired, if error thrown here file truncated
                    sam_line2 = sam.readline()
                    sam_dict2 = self.process_sam_line(sam_line2)
                    if sam_id != sam_dict1['QNAME']:
                        sam_id = sam_dict1['QNAME']
                        yield {f'{sam_dict1["QNAME"]}_1': sam_dict1, f'{sam_dict2["QNAME"]}_2': sam_dict2}
                else:
                    if sam_id != sam_dict1['QNAME']:
                        sam_id = sam_dict1['QNAME']
                        yield {sam_dict1['QNAME']: sam_dict1}

    @staticmethod
    def process_sam_line(sam_line):
        if sam_line:
            processed_sam_line: list = sam_line.replace('\n', '').split('\t')
            QNAME, original_sequence = processed_sam_line[0].split('_BSBolt_')
            FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = processed_sam_line[1:11]
            SAM_TAGS = None
            if len(sam_line) > 11:
                SAM_TAGS = processed_sam_line[11:]
            return dict(QNAME=QNAME, FLAG=FLAG, RNAME=RNAME, POS=POS, MAPQ=MAPQ, CIGAR=CIGAR, RNEXT=RNEXT, PNEXT=PNEXT,
                        TLEN=TLEN, SEQ=SEQ, QUAL=QUAL, SAM_TAGS=SAM_TAGS, original_sequence=original_sequence)
