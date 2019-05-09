from BSB.BSB_Utils.FastaIterator import OpenFasta



def get_reference_contigs(reference_file):
    """Retrieve contig sequence from reference fasta rather than index to ensure simulation will work without
    a generated index. """
    # store contig_id and matching sequencing
    contig_dict = {}
    contig_id = None
    contig_sequence = []
    for contig_label, line in OpenFasta(reference_file):
        if contig_label:
            # join sequence and add to dict if new contig encounted
            if contig_id:
                contig_dict[contig_id] = ''.join(contig_sequence)
                contig_sequence = []
            contig_id = line.replace('>', '').split()[0]
        else:
            contig_sequence.append(line)
    # join contig and add to dict after iteration finishes
    contig_dict[contig_id] = ''.join(contig_sequence)
    return contig_dict

