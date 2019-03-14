import networkx as nx


def build_alignment_graph(assembly_graph, alignment_blocks):
    alignment_graph = nx.MultiDiGraph()
    node_curr_id = max(assembly_graph.nodes) + 1
    distances_between_blocks = dict()

    for (node_from, node_to, data) in assembly_graph.edges(data=True):
        edge_blocks = alignment_blocks.get(data["name"])

        if edge_blocks is None:
            if node_from != node_to:
                alignment_graph.add_edge(node_from, node_to, **data)
            continue

        if len(edge_blocks) == 1:
            alignment_graph.add_edge(node_from, node_to, name=edge_blocks[0].signed_id())
            continue

        edge_blocks.sort(key=lambda block: block.start)
        alignment_graph.add_edge(node_from, node_curr_id, name=edge_blocks[0].signed_id())
        node_curr_id += 1

        for i in range(1, len(edge_blocks) - 1):
            alignment_graph.add_edge(node_curr_id - 1, node_curr_id, name=edge_blocks[i].signed_id())

            dist = edge_blocks[i].start - edge_blocks[i - 1].end
            distances_between_blocks[node_curr_id - 1] = {"distance": dist}

            node_curr_id += 1

        alignment_graph.add_edge(node_curr_id - 1, node_to, name=edge_blocks[-1].signed_id())

        dist = edge_blocks[-1].start - edge_blocks[-2].end
        distances_between_blocks[node_curr_id - 1] = {"distance": dist}

        node_curr_id += 1

        if node_from == node_to:
            dist = edge_blocks[-1].seq_length - edge_blocks[-1].end + edge_blocks[0].start
            distances_between_blocks[node_to] = {"distance": dist}

    nx.set_node_attributes(alignment_graph, distances_between_blocks)
    return alignment_graph
