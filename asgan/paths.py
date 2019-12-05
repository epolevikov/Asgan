import networkx as nx
from src.synteny_blocks import SequenceBlock
from src.adjacency_graph import build_contracted_adjacency_graph
from src.common import build_id2block_dict, build_block2edge_dict


def build_synteny_paths(path_components):
    components = nx.connected_component_subgraphs(path_components)
    node_labels = nx.get_node_attributes(path_components, "label")

    paths = []
    for component in components:
        component_is_processed = False
        path = None

        for node in component.nodes():
            if component.degree[node] == 1:
                path = list(nx.dfs_preorder_nodes(component, node))
                component_is_processed = True
                break

        if not component_is_processed:
            for (node_from, node_to) in component.edges():
                node_from_label = node_labels[node_from]
                node_to_label = node_labels[node_to]

                if node_from_label[:-1] != node_to_label[:-1]:
                    component.remove_edge(node_from, node_to)
                    path = list(nx.dfs_preorder_nodes(component, node_from))
                    # path.append(node_from)
                    break

        paths.append(path)

    synteny_paths = []
    for path in paths:
        synteny_path = []

        for i in range(0, len(path), 2):
            node_label = node_labels[path[i]]
            block_id = ["+", "-"][node_label[-1] == "h"] + node_label[:-1]
            synteny_path.append(block_id)

        synteny_paths.append(synteny_path)

    return synteny_paths


def build_path_sequences(synteny_blocks, synteny_paths, adjacency_graph):
    contracted_adjacency_graph = build_contracted_adjacency_graph(adjacency_graph)
    edge2data = build_edge2data_dict(contracted_adjacency_graph)
    block2edge = build_block2edge_dict(adjacency_graph)
    id2block = build_id2block_dict(synteny_blocks)

    path_sequences = []
    for synteny_path in synteny_paths:
        path_sequence = build_path_sequence(synteny_path, contracted_adjacency_graph,
                                            id2block, block2edge, edge2data)
        path_sequences.append(path_sequence)

    return path_sequences


def build_path_sequence(synteny_path, contracted_adjacency_graph,
                        id2block, block2edge, edge2data):
    if len(synteny_path) == 1:
        return [id2block[synteny_path[0]]]

    path_sequence = []

    block_from = id2block[synteny_path[0]]
    block_to = id2block[synteny_path[1]]
    path_between_blocks = build_path_between_blocks(block_from, block_to,
                                                    contracted_adjacency_graph,
                                                    block2edge, edge2data)

    path_sequence.append(block_from)
    path_sequence.append(path_between_blocks)
    path_sequence.append(block_to)

    for i in range(2, len(synteny_path)):
        block_from = block_to
        block_to = id2block[synteny_path[i]]
        path_between_blocks = build_path_between_blocks(block_from, block_to,
                                                        contracted_adjacency_graph,
                                                        block2edge, edge2data)

        path_sequence.append(path_between_blocks)
        path_sequence.append(block_to)

    return path_sequence


def build_path_between_blocks(block_from, block_to,
                              contracted_adjacency_graph,
                              block2edge, edge2data):
    (from_start, from_end) = block2edge[block_from.signed_id()]
    (to_start, to_end) = block2edge[block_to.signed_id()]

    if from_end == to_start:
        dist = contracted_adjacency_graph.nodes[from_end].get("distance")

        if dist is None or block_from.sequence_name != block_to.sequence_name:
            block1 = SequenceBlock(id=None, sequence_name=block_from.sequence_name,
                                   sequence_length=block_from.sequence_length,
                                   start=block_from.end, end=block_from.sequence_length)

            block2 = SequenceBlock(id=None, sequence_name=block_to.sequence_name,
                                   sequence_length=block_to.sequence_length,
                                   start=0, end=block_to.start)

            return [block1, block2]
        else:
            block = SequenceBlock(id=None, sequence_name=block_from.sequence_name,
                                  sequence_length=block_from.sequence_length,
                                  start=block_from.end,
                                  end=block_to.start
                                  if block_to.start > block_from.start
                                  else block_to.sequence_length)  # in case if the sequence is circular

            return [block]

    path_nodes = nx.dijkstra_path(contracted_adjacency_graph, from_end, to_start)

    path = []
    if not path_nodes:
        data = edge2data[(from_end, to_start)]
        block = SequenceBlock(id=None, sequence_name=data["name"],
                              sequence_length=data["weight"],
                              start=0, end=data["weight"])
        path.append(block)
        return path

    for i in range(len(path_nodes) - 1):
        node_from, node_to = path_nodes[i], path_nodes[i + 1]
        data = edge2data[(node_from, node_to)]
        block = SequenceBlock(id=None, sequence_name=data["name"],
                              sequence_length=data["weight"],
                              start=0, end=data["weight"])
        path.append(block)

    return path


def build_edge2data_dict(contracted_adjacency_graph):
    edge2data_dict = dict()

    for (node_from, node_to, data) in contracted_adjacency_graph.edges(data=True):
        prev_data = edge2data_dict.get((node_from, node_to))

        if prev_data is None or prev_data["weight"] > data["weight"]:
            edge2data_dict[(node_from, node_to)] = data

    return edge2data_dict
