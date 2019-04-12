from asgan.output_generator import pretty_number


class SyntenyBlock:
    def __init__(self, id, sequence_name, sequence_length, start, end):
        self.id = id
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
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


def extract_synteny_blocks(hits):
    synteny_blocks_query, synteny_blocks_target = [], []

    for hit in hits:
        synteny_block_query = SyntenyBlock(hit.id,
                                           hit.query_name, hit.query_len,
                                           hit.query_start, hit.query_end)
        synteny_block_target = SyntenyBlock(hit.id,
                                            hit.target_name, hit.target_len,
                                            hit.target_start, hit.target_end)

        synteny_blocks_query.append(synteny_block_query)
        synteny_blocks_target.append(synteny_block_target)

    return synteny_blocks_query, synteny_blocks_target


def group_by_sequence(synteny_blocks):
    grouped_synteny_blocks = dict()

    for synteny_block in synteny_blocks:
        if synteny_block.sequence_name not in grouped_synteny_blocks:
            grouped_synteny_blocks[synteny_block.sequence_name] = []

        grouped_synteny_blocks[synteny_block.sequence_name].append(synteny_block)

    return grouped_synteny_blocks


def set_block_attributes(paths):
    colors = ["green", "blue", "gold", "red", "purple", "darkorange",
              "hotpink", "khaki", "lightblue", "thistle", "tan"]
    block_colors = dict()
    block_styles = dict()
    i = 0

    paths.sort(key=lambda path: len(path), reverse=True)

    while i < len(colors) and i < len(paths):
        for block in paths[i]:
            block_colors[block] = colors[i]
            block_styles[block] = "solid"

            block_colors[inv_block(block)] = colors[i]
            block_styles[inv_block(block)] = "dashed"

        i += 1

    while i < len(paths):
        for block in paths[i]:
            block_colors[block] = "gray"
            block_styles[block] = "solid"

            block_colors[inv_block(block)] = "gray"
            block_styles[inv_block(block)] = "dashed"

        i += 1

    return {"color": block_colors, "style": block_styles}


def inv_block(block_id):
    return ["+", "-"][block_id[0] == "+"] + block_id[1:]


'''
def build_from_sequences(assembly_graph, repeats):
    sequences = []

    for (_, _, data) in assembly_graph.edges(data=True):
        if data["name"][:-1] not in repeats:
            sequences.append((data["name"], data["length"]))

    sequences.sort(key=lambda p: p[0])
    alignment_blocks = dict()

    block_id = 1
    for i, (sequence_name, sequence_length) in enumerate(sequences):
        block = AlignmentBlock(block_id, sequence_name, sequence_length,
                               0, sequence_length)
        alignment_blocks[sequence_name] = [block]

        if (i + 1) % 2 == 1:
            block_id = -block_id
        else:
            block_id = -block_id + 1

    return alignment_blocks
'''
