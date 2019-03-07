import networkx as nx
from asgan.utils import DisjointSet


class Sequence:
    def __init__(self, name, length, strand):
        self.name = name + strand
        self.length = length


'''
class Edge:
    def __init__(self, name, length, node_from, node_to):
        self.name = name
        self.length = length
        self.node_from = node_from
        self.node_to = node_to

    def __str__(self):
        return "name={}, length={}, node_from={}, node_to={}".format(
            self.name, self.length, self.node_from, self.node_to)
'''


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
    assembly_graph.add_nodes_from(disjoint_set.get_unique_parents())

    for i, seq in enumerate(sequences):
        node_from = disjoint_set.find(2 * i)
        node_to = disjoint_set.find(2 * i + 1)
        assembly_graph.add_edge(node_from, node_to,
                                name=seq.name,
                                length=seq.length)

    return assembly_graph

    '''
    assembly_graph = nx.MultiDiGraph()

    for v in range(len(digraph)):
        is_self_loop = (v in digraph[v])
        assembly_graph.add_edge(2 * v, 2 * v + 1,
                                name="{}".format(sequences[v].name),
                                self_loop=is_self_loop)

    for v, neighbour in enumerate(digraph):
        for u in neighbour:
            if v != u:
                assembly_graph.add_edge(2 * v + 1, 2 * u,
                                        name="dummy",
                                        key="dummy")

    has_dummy_edges = True
    while has_dummy_edges:
        has_dummy_edges = False
        for (node_from, node_to, data) in assembly_graph.edges(data=True):
            if data["name"] == "dummy":
                if assembly_graph.number_of_edges(node_from, node_to) > 1:
                    assembly_graph.remove_edge(node_from, node_to,
                                               key="dummy")
                else:
                    assembly_graph = nx.contracted_edge(assembly_graph,
                                                        (node_from, node_to),
                                                        self_loops=False)
                has_dummy_edges = True
                break

    has_self_loops = True
    while has_self_loops:
        has_self_loops = False
        for (node_from, node_to, data) in assembly_graph.edges(data=True):
            if node_from != node_to and data["self_loop"]:
                assembly_graph = nx.contracted_edge(assembly_graph,
                                                    (node_from, node_to))
                has_self_loops = True
                break

    edges = build_edge_graph(digraph, sequences)

    assembly_graph = nx.MultiDiGraph()
    for edge in edges:
        assembly_graph.add_edge(edge.node_from, edge.node_to,
                                name=edge.name, length=edge.length)

    return assembly_graph
    '''


'''
def build_edge_graph(digraph, sequences):
    def unite_nodes(i, j):
        nodes_union = nodes[i].union(nodes[j])
        for node in nodes_union:
            nodes[node] = nodes_union

    edges = []
    node_ctr = 0

    for node_id in range(len(digraph)):
        name, length = sequences[node_id].name, sequences[node_id].length
        node_from, node_to = node_ctr, node_ctr + 1
        edges.append(Edge(name, length, node_from, node_to))
        node_ctr += 2

    nodes = {node_id: {node_id} for node_id in range(2 * len(edges))}

    for node, neighbour in enumerate(digraph):
        for u in neighbour:
            unite_nodes(edges[node].node_to, edges[u].node_from)

    nodes = {node_id: min(nodes[node_id]) for node_id in nodes}

    for edge in edges:
        edge.node_from = nodes[edge.node_from]
        edge.node_to = nodes[edge.node_to]

    return edges
'''


def remove_components_without_alignment_blocks(assembly_graph,
                                               alignment_blocks):
    def check_component(wcc):
        for (node_from, node_to, data) in wcc.edges(data=True):
            if data["name"] in alignment_blocks:
                return True

        return False

    components_with_alignment_blocks = nx.MultiDiGraph()
    for wcc in nx.weakly_connected_component_subgraphs(assembly_graph,
                                                       copy=False):
        if check_component(wcc):
            components_with_alignment_blocks = \
                nx.union(components_with_alignment_blocks, wcc)

    return components_with_alignment_blocks
