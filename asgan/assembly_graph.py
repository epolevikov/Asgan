import networkx as nx
from asgan.utils import DisjointSet


class Sequence:
    def __init__(self, name, length, strand, depth):
        self.name = name + strand
        self.length = length
        self.depth = depth


def build_assembly_graph(raw_sequences, links):
    sequences = list()
    sequence2id = dict()
    curr_id = 0

    for seq in raw_sequences:
        sequences.append(Sequence(seq.name, seq.length, "+", seq.depth))
        sequence2id[seq.name + "+"] = curr_id

        sequences.append(Sequence(seq.name, seq.length, "-", seq.depth))
        sequence2id[seq.name + "-"] = curr_id + 1

        curr_id += 2

    disjoint_set = DisjointSet(2 * len(sequences))

    for link in links:
        id_from = sequence2id[link.from_name + link.from_strand]
        id_to = sequence2id[link.to_name + link.to_strand]
        disjoint_set.union(2 * id_from + 1, 2 * id_to)

    assembly_graph = nx.MultiDiGraph()

    for i, seq in enumerate(sequences):
        node_from = disjoint_set.find(2 * i)
        node_to = disjoint_set.find(2 * i + 1)
        assembly_graph.add_edge(node_from, node_to, name=seq.name,
                                length=seq.length, depth=seq.depth,
                                is_repeat=False)

    return assembly_graph


'''
def mark_repeats(assembly_graph, normalize_depth=False):
    if not normalize_depth:
        for (_, _, data) in assembly_graph.edges(data=True):
            is_repeat = (data["depth"] > 1 or data["length"] < 50000)
            data["is_repeat"] = is_repeat
    else:
        depths, lengths = [], []

        for (_, _, data) in assembly_graph.edges(data=True):
            depths.append(data["depth"])
            lengths.append(data["length"])

        weighted_depths_sum = sum([depths[i] * lengths[i] for i in range(len(depths))])
        weighted_mean_depth = weighted_depths_sum / sum(lengths)

        for (_, _, data) in assembly_graph.edges(data=True):
            is_repeat = (data["depth"] > 2 * weighted_mean_depth or data["length"] < 50000)
            data["is_repeat"] = is_repeat


def mark_sequences_without_alignment_blocks(assembly_graph, alignment_blocks):
    for (_, _, data) in assembly_graph.edges(data=True):
        if data["name"] not in alignment_blocks:
            data["is_repeat"] = True


def get_repeats(assembly_graph):
    repeats = set()

    for (_, _, data) in assembly_graph.edges(data=True):
        if data["is_repeat"]:
            repeats.add(data["name"][:-1])

    return repeats
'''
