import networkx as nx


def build_adjacency_graph(assembly_graph, synteny_blocks):
    adjacency_graph = nx.MultiDiGraph()
    node_curr_id = max(assembly_graph.nodes) + 1
    distances_between_blocks = dict()

    for (node_from, node_to, data) in assembly_graph.edges(data=True):
        edge_synteny_blocks = synteny_blocks.get(data["name"])

        if edge_synteny_blocks is None:
            if node_from != node_to:
                adjacency_graph.add_edge(node_from, node_to, **data)
            continue

        if len(edge_synteny_blocks) == 1:
            edge_name = edge_synteny_blocks[0].signed_id()
            adjacency_graph.add_edge(node_from, node_to, name=edge_name)
            continue

        edge_synteny_blocks.sort(key=lambda block: block.start)

        edge_name = edge_synteny_blocks[0].signed_id()
        adjacency_graph.add_edge(node_from, node_curr_id, name=edge_name)
        node_curr_id += 1

        for i in range(1, len(edge_synteny_blocks) - 1):
            edge_name = edge_synteny_blocks[i].signed_id()
            adjacency_graph.add_edge(node_curr_id - 1, node_curr_id, name=edge_name)

            dist = edge_synteny_blocks[i].start - edge_synteny_blocks[i - 1].end
            distances_between_blocks[node_curr_id - 1] = {"distance": dist}

            node_curr_id += 1

        edge_name = edge_synteny_blocks[-1].signed_id()
        adjacency_graph.add_edge(node_curr_id - 1, node_to, name=edge_name)

        dist = edge_synteny_blocks[-1].start - edge_synteny_blocks[-2].end
        distances_between_blocks[node_curr_id - 1] = {"distance": dist}

        node_curr_id += 1

        if node_from == node_to:
            dist = edge_synteny_blocks[-1].sequence_length - edge_synteny_blocks[-1].end
            dist += edge_synteny_blocks[0].start
            distances_between_blocks[node_to] = {"distance": dist}

    nx.set_node_attributes(adjacency_graph, distances_between_blocks)
    return adjacency_graph
