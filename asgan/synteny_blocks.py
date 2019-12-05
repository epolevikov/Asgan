from src.output_generator import pretty_number


class SequenceBlock:
    def __init__(self, id, sequence_name, sequence_length, start, end):
        self.id = id
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
        self.start = start
        self.end = end

    def signed_id(self):
        if self.id is None:
            return ""

        return ["-", "+"][self.id > 0] + str(abs(self.id))

    def length(self):
        return self.end - self.start

    def __str__(self):
        name = self.sequence_name
        start = pretty_number(self.start)
        end = pretty_number(self.end)
        id = "-"

        if self.id is not None:
            id = self.signed_id()

        return "{}\t{}\t{}\t{}".format(name, start, end, id)


class SyntenyBlock(SequenceBlock):
    pass


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

    return (group_by_sequence(synteny_blocks_query),
            group_by_sequence(synteny_blocks_target))


def group_by_sequence(synteny_blocks):
    grouped_synteny_blocks = dict()

    for synteny_block in synteny_blocks:
        if synteny_block.sequence_name not in grouped_synteny_blocks:
            grouped_synteny_blocks[synteny_block.sequence_name] = []

        grouped_synteny_blocks[synteny_block.sequence_name].append(synteny_block)

    return grouped_synteny_blocks


def set_block_attributes(paths):
    colors = ["blue", "green", "gold", "red", "purple", "darkorange",
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
