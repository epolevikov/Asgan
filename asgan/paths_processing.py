import networkx as nx


def reconstruct_paths(breakpoint_graph, max_matching):
    paths = nx.Graph()

    for node, data in breakpoint_graph.nodes(data=True):
        paths.add_node(node, **data)

    for (node_from, node_to) in max_matching:
        paths.add_edge(node_from, node_to)

    for node_from, data_from in breakpoint_graph.nodes(data=True):
        for node_to, data_to in breakpoint_graph.nodes(data=True):
            if node_from == node_to:
                continue

            label_from, label_to = data_from["label"], data_to["label"]
            if label_from[:-1] == label_to[:-1]:
                paths.add_edge(node_from, node_to)

    return paths
