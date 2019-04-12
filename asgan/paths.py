import networkx as nx


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
                    path.append(node_from)
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

    path_sequence.extend([block_from, path_between_blocks, block_to])

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

        if dist is None:
            dist = block_to.start + (block_from.seq_length - block_from.end)

        return [("", dist)]

    path_nodes = nx.dijkstra_path(contracted_adjacency_graph, from_end, to_start)

    path = []
    if not path_nodes:
        data = edge2data[(from_end, to_start)]
        path.append((data["name"], data["weight"]))
        return path

    for i in range(len(path_nodes) - 1):
        node_from, node_to = path_nodes[i], path_nodes[i + 1]
        data = edge2data[(node_from, node_to)]
        path.append((data["name"], data["weight"]))

    return path


def build_contracted_adjacency_graph(adjacency_graph):
    contracted_adjacency_graph = nx.MultiDiGraph()

    for (node, data) in adjacency_graph.nodes(data=True):
        contracted_adjacency_graph.add_node(node, **data)

    for (node_from, node_to, data) in adjacency_graph.edges(data=True):
        if data["name"].startswith("+") or data["name"].startswith("-"):
            continue

        contracted_adjacency_graph.add_edge(node_from, node_to,
                                            name=data["name"],
                                            weight=data["length"])

    return contracted_adjacency_graph


def build_edge2data_dict(contracted_adjacency_graph):
    edge2data_dict = dict()

    for (node_from, node_to, data) in contracted_adjacency_graph.edges(data=True):
        prev_data = edge2data_dict.get((node_from, node_to))

        if prev_data is None or prev_data["weight"] > data["weight"]:
            edge2data_dict[(node_from, node_to)] = data

    return edge2data_dict


def build_block2edge_dict(adjacency_graph):
    block2edge = dict()

    for (node_from, node_to, data) in adjacency_graph.edges(data=True):
        edge_name = data["name"]

        if edge_name.startswith("+") or edge_name.startswith("-"):
            block2edge[edge_name] = (node_from, node_to)

    return block2edge


def build_id2block_dict(synteny_blocks):
    id2block = dict()

    for blocks in synteny_blocks.values():
        for block in blocks:
            id2block[block.signed_id()] = block

    return id2block


'''
from asgan.breakpoint_graph import build_block2edge_dict, build_contracted_alignment_graph


def reconstruct_alignment_block_paths(alignment_graph_query, alignment_blocks_query,
                                      alignment_graph_target, alignment_blocks_target,
                                      path_components):
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
                    path.append(node_from)
                    break

        paths.append(path)

    alignment_block_paths = []
    for path in paths:
        alignment_block_path = []

        for i in range(0, len(path), 2):
            node_label = node_labels[path[i]]
            alignment_block = ["+", "-"][node_label[-1] == "h"] + node_label[:-1]
            alignment_block_path.append(alignment_block)

        complement_alignment_block_path = _complement_path(alignment_block_path)
        alignment_block_paths.append((alignment_block_path, complement_alignment_block_path))

    return alignment_block_paths


def reconstruct_full_paths(alignment_block_paths, alignment_graph, alignment_blocks):
    contracted_alignment_graph = build_contracted_alignment_graph(alignment_graph)
    edge2data = _build_edge2data_dict(contracted_alignment_graph)
    block2edge = build_block2edge_dict(alignment_graph)
    id2block = build_id2block_dict(alignment_blocks)

    full_paths = []
    for alignment_block_path in alignment_block_paths:
        full_path = reconstruct_full_path_from_alignment_block_path(
            alignment_block_path[0], contracted_alignment_graph,
            id2block, block2edge, edge2data)

        full_complement_path = reconstruct_full_path_from_alignment_block_path(
            alignment_block_path[1], contracted_alignment_graph,
            id2block, block2edge, edge2data)

        full_paths.append((full_path, full_complement_path))

    return full_paths


def reconstruct_full_path_from_alignment_block_path(alignment_block_path,
                                                    contracted_alignment_graph,
                                                    id2block, block2edge, edge2data):
    if len(alignment_block_path) == 1:
        return [id2block[alignment_block_path[0]]]

    full_path = []

    block_from = id2block[alignment_block_path[0]]
    block_to = id2block[alignment_block_path[1]]
    path_between_blocks = reconstruct_path_between_blocks(
        block_from, block_to, contracted_alignment_graph,
        block2edge, edge2data)

    full_path.extend([block_from, path_between_blocks, block_to])

    for i in range(2, len(alignment_block_path)):
        block_from = block_to
        block_to = id2block[alignment_block_path[i]]
        path_between_blocks = reconstruct_path_between_blocks(
            block_from, block_to, contracted_alignment_graph,
            block2edge, edge2data)

        full_path.append(path_between_blocks)
        full_path.append(block_to)

    return full_path


def reconstruct_path_between_blocks(block_from, block_to, contracted_alignment_graph, block2edge, edge2data):
    (from_start, from_end) = block2edge[block_from.signed_id()]
    (to_start, to_end) = block2edge[block_to.signed_id()]

    if from_end == to_start:
        dist = contracted_alignment_graph.nodes[from_end].get("distance")

        if dist is None:
            dist = block_to.start + (block_from.seq_length - block_from.end)

        return [("", dist)]

    path_nodes = nx.dijkstra_path(contracted_alignment_graph, from_end, to_start)

    path = []
    if not path_nodes:
        data = edge2data[(from_end, to_start)]
        path.append((data["name"], data["weight"]))
        return path

    for i in range(len(path_nodes) - 1):
        node_from, node_to = path_nodes[i], path_nodes[i + 1]
        data = edge2data[(node_from, node_to)]
        path.append((data["name"], data["weight"]))

    return path


def build_id2block_dict(alignment_blocks):
    id2block = dict()

    for blocks in alignment_blocks.values():
        for block in blocks:
            id2block[block.signed_id()] = block

    return id2block


def _build_edge2data_dict(contracted_alignment_graph):
    edge2data_dict = dict()

    for (node_from, node_to, data) in contracted_alignment_graph.edges(data=True):
        prev_data = edge2data_dict.get((node_from, node_to))

        if prev_data is None or prev_data["weight"] > data["weight"]:
            edge2data_dict[(node_from, node_to)] = data

    return edge2data_dict


def _complement_path(path):
    comp_path = []

    for label in reversed(path):
        comp_label = ["+", "-"][label[0] == "+"] + label[1:]
        comp_path.append(comp_label)

    return comp_path
'''
