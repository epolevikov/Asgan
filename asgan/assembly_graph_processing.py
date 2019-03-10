import networkx as nx
from asgan.utils import DisjointSet


class Sequence:
    def __init__(self, name, length, strand):
        self.name = name + strand
        self.length = length


def build_assembly_graph(raw_sequences, links):
    sequences = list()
    sequence2id = dict()
    curr_id = 0

    for seq in raw_sequences:
        sequences.append(Sequence(seq.name, seq.length, "+"))
        sequence2id[seq.name + "+"] = curr_id

        sequences.append(Sequence(seq.name, seq.length, "-"))
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
        assembly_graph.add_edge(node_from, node_to, name=seq.name, length=seq.length)

    return assembly_graph


def remove_components_without_alignment_blocks(assembly_graph, alignment_blocks):
    def check_component(wcc):
        for (node_from, node_to, data) in wcc.edges(data=True):
            if data["name"] in alignment_blocks:
                return True

        return False

    components_with_alignment_blocks = nx.MultiDiGraph()
    for wcc in nx.weakly_connected_component_subgraphs(assembly_graph, copy=False):
        if check_component(wcc):
            components_with_alignment_blocks = nx.union(components_with_alignment_blocks, wcc)

    return components_with_alignment_blocks
