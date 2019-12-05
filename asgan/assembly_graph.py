import networkx as nx
from asgan.utils import DisjointSet
from asgan.gfa_parser import parse_gfa


class Sequence:
    def __init__(self, name, length, strand, is_repeat):
        self.name = name + strand
        self.length = length
        self.is_repeat = is_repeat


def parse_assembly_graph(gfa_file):
    sequences, links = parse_gfa(gfa_file)
    return build(sequences, links)


def build(raw_sequences, links):
    sequences = list()
    sequence2id = dict()
    curr_id = 0

    for seq in raw_sequences:
        sequences.append(Sequence(seq.name, seq.length, "+", seq.is_repeat))
        sequence2id[seq.name + "+"] = curr_id

        sequences.append(Sequence(seq.name, seq.length, "-", seq.is_repeat))
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
                                length=seq.length,
                                is_repeat=seq.is_repeat)

    return assembly_graph


def get_repeats(assembly_graph):
    repeats = set()

    for (_, _, data) in assembly_graph.edges(data=True):
        if data["is_repeat"]:
            repeats.add(data["name"][:-1])

    return repeats
