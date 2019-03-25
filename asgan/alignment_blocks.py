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

    def length(self):
        return self.end - self.start + 1

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


def set_block_attributes(paths):
    colors = ["green", "blue", "gold", "red", "purple", "darkorange",
              "hotpink", "khaki", "lightblue", "thistle", "tan"]
    block_colors = dict()
    block_styles = dict()
    i = 0

    paths.sort(key=lambda pair: len(pair[0]), reverse=True)

    while i < len(colors) and i < len(paths):
        for block in paths[i][0]:
            block_colors[block] = colors[i]
            block_styles[block] = "solid"

        for block in paths[i][1]:
            block_colors[block] = colors[i]
            block_styles[block] = "dashed"

        i += 1

    while i < len(paths):
        for block in paths[i][0]:
            block_colors[block] = "gray"
            block_styles[block] = "solid"

        for block in paths[i][1]:
            block_colors[block] = "gray"
            block_styles[block] = "dashed"

        i += 1

    return block_colors, block_styles
