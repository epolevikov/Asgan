import networkx as nx


def reconstruct_paths(breakpoint_graph, max_matching):
    paths_graph = _build_paths_graph(breakpoint_graph, max_matching)
    paths = []
    components = nx.connected_component_subgraphs(paths_graph, copy=False)

    for component in components:
        component_in_processed = False
        path = None

        for node in component.nodes():
            if component.degree[node] == 1:
                path = nx.dfs_preorder_nodes(component, node)
                component_in_processed = True
                break

        if not component_in_processed:
            path = nx.dfs_preorder_nodes(component)

        paths.append(list(path))

    node_labels = nx.get_node_attributes(paths_graph, "label")
    labeled_paths = []

    for path in paths:
        labeled_path = []

        for i in range(0, len(path), 2):
            node_label = node_labels[path[i]]
            node_label = ["+", "-"][node_label[-1] == "h"] + node_label[:-1]
            labeled_path.append(node_label)

        complement_labeled_path = _complement_path(labeled_path)
        labeled_paths.append([labeled_path, complement_labeled_path])

    return paths_graph, labeled_paths


def _build_paths_graph(breakpoint_graph, max_matching):
    paths_graph = nx.Graph()

    for node, data in breakpoint_graph.nodes(data=True):
        paths_graph.add_node(node, **data)

    for (node_from, node_to) in max_matching:
        paths_graph.add_edge(node_from, node_to)

    for node_from, data_from in breakpoint_graph.nodes(data=True):
        for node_to, data_to in breakpoint_graph.nodes(data=True):
            if node_from == node_to:
                continue

            label_from, label_to = data_from["label"], data_to["label"]
            if label_from[:-1] == label_to[:-1]:
                paths_graph.add_edge(node_from, node_to)

    return paths_graph


def _complement_path(path):
    comp_path = []

    for label in reversed(path):
        comp_label = ["+", "-"][label[0] == "+"] + label[1:]
        comp_path.append(comp_label)

    return comp_path
