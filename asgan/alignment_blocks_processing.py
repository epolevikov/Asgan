from asgan.output_generator import pretty_number


class AlignmentBlock:
    def __init__(self, id, seq_name, seq_length, start, end):
        self.id = id
        self.seq_name = seq_name
        self.seq_length = seq_length
        self.start = start
        self.end = end

    def signed_id(self):
        return ["-", "+"][self.id > 0] + str(abs(self.id))

    def start_gap_length(self):
        return self.start

    def end_gap_length(self):
        return self.seq_length - self.end

    def __str__(self):
        return "{}\t{}\t{}".format(self.signed_id(),
                                   pretty_number(self.start),
                                   pretty_number(self.end))


def extract_alignment_blocks(hits):
    aln_blocks_query, aln_blocks_target = [], []

    for hit in hits:
        aln_block_query = AlignmentBlock(hit.id,
                                         hit.query_name, hit.query_len,
                                         hit.query_start, hit.query_end)
        aln_block_target = AlignmentBlock(hit.id,
                                          hit.target_name, hit.target_len,
                                          hit.target_start, hit.target_end)
        aln_blocks_query.append(aln_block_query)
        aln_blocks_target.append(aln_block_target)

    return aln_blocks_query, aln_blocks_target


def group_by_sequence(alignment_blocks):
    grouped_aln_blocks = dict()

    for alignment_block in alignment_blocks:
        if alignment_block.seq_name not in grouped_aln_blocks:
            grouped_aln_blocks[alignment_block.seq_name] = []

        grouped_aln_blocks[alignment_block.seq_name].append(alignment_block)

    return grouped_aln_blocks
