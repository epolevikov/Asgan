import networkx as nx


class Sequence:
    def __init__(self, name, length, strand):
        self.name = name + strand
        self.length = length


class Edge:
    def __init__(self, name, length, node_from, node_to):
        self.name = name
        self.length = length
        self.node_from = node_from
        self.node_to = node_to

    def __str__(self):
        return "name={}, length={}, node_from={}, node_to={}".format(
            self.name, self.length, self.node_from, self.node_to)


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

    n_vertices = len(sequences)
    digraph = [[] for _ in range(n_vertices)]

    for link in links:
        id_from = sequence2id[link.from_name + link.from_strand]
        id_to = sequence2id[link.to_name + link.to_strand]
        digraph[id_from].append(id_to)

    edges = build_edge_graph(digraph, sequences)

    assembly_graph = nx.MultiDiGraph()
    for edge in edges:
        assembly_graph.add_edge(edge.node_from, edge.node_to,
                                name=edge.name, length=edge.length)

    return assembly_graph


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


def remove_components_without_alignment_blocks(assembly_graph,
                                               alignment_blocks):
    def check_component(wcc):
        for (node_from, node_to, data) in wcc.edges(data=True):
            if data["name"] in alignment_blocks:
                return True

        return False

    components_with_alignment_blocks = nx.MultiDiGraph()
    for wcc in nx.weakly_connected_component_subgraphs(assembly_graph):
        if check_component(wcc):
            components_with_alignment_blocks = \
                nx.union(components_with_alignment_blocks, wcc)

    return components_with_alignment_blocks
