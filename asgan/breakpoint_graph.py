import networkx as nx


def build_breakpoint_graph(adjacency_graph_query, synteny_blocks_query,
                           adjacency_graph_target, synteny_blocks_target):
    number_synteny_blocks = count_number_synteny_blocks(adjacency_graph_query)
    breakpoint_graph = nx.Graph()
    labels = dict()

    for i in range(1, number_synteny_blocks + 1):
        tail, head = "{}t".format(i), "{}h".format(i)

        breakpoint_graph.add_node(i, label=tail)
        breakpoint_graph.add_node(number_synteny_blocks + i, label=head)

        labels[tail] = i
        labels[head] = number_synteny_blocks + i

    contracted_adjacency_graph_query = build_contracted_adjacency_graph(adjacency_graph_query)
    contracted_adjacency_graph_target = build_contracted_adjacency_graph(adjacency_graph_target)

    block2edge_query = build_block2edge_dict(adjacency_graph_query)
    block2edge_target = build_block2edge_dict(adjacency_graph_target)

    id2block_query = build_id2block_dict(synteny_blocks_query)
    id2block_target = build_id2block_dict(synteny_blocks_target)

    signs = [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]

    for i in range(1, number_synteny_blocks):
        for j in range(i + 1, number_synteny_blocks + 1):
            for (sign1, sign2) in signs:
                block_fwd_from, block_fwd_to = make_signed_blocks(i, j, sign1, sign2)
                block_inv_from, block_inv_to = make_signed_blocks(j, i, inv_sign(sign2),
                                                                  inv_sign(sign1))

                if not check_adjacency(block_fwd_from, block_fwd_to,
                                       contracted_adjacency_graph_query,
                                       block2edge_query, id2block_query):
                    continue

                if not check_adjacency(block_inv_from, block_inv_to,
                                       contracted_adjacency_graph_query,
                                       block2edge_query, id2block_query):
                    continue

                if not check_adjacency(block_fwd_from, block_fwd_to,
                                       contracted_adjacency_graph_target,
                                       block2edge_target, id2block_target):
                    continue

                if not check_adjacency(block_inv_from, block_inv_to,
                                       contracted_adjacency_graph_target,
                                       block2edge_target, id2block_target):
                    continue

                label_from = "{}{}".format(i, ["t", "h"][sign1 == "+"])
                label_to = "{}{}".format(j, ["h", "t"][sign2 == "+"])
                node_from = labels[label_from]
                node_to = labels[label_to]
                breakpoint_graph.add_edge(node_from, node_to)

    return breakpoint_graph


def count_number_synteny_blocks(adjacency_graph):
    number_synteny_blocks = 0

    for (node_from, node_to, data) in adjacency_graph.edges(data=True):
        if data["name"].startswith("+") or data["name"].startswith("-"):
            number_synteny_blocks += 1

    return number_synteny_blocks // 2


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


def make_signed_blocks(id_from, id_to, sign_from, sign_to):
    signed_id_from = "{}{}".format(sign_from, id_from)
    signed_id_to = "{}{}".format(sign_to, id_to)
    return signed_id_from, signed_id_to


def inv_sign(sign):
    return ["+", "-"][sign == "+"]


def check_adjacency(id_from, id_to, contracted_adjacency_graph, block2edge, id2block):
    from_start, from_end = block2edge[id_from]
    to_start, to_end = block2edge[id_to]

    if from_end == to_start:
        dist = contracted_adjacency_graph.nodes[from_end].get("distance")
        block_from, block_to = id2block[id_from], id2block[id_to]

        if dist is None or block_from.sequence_name != block_to.sequence_name:
            dist = block_to.start + (block_from.sequence_length - block_from.end)

        return dist < 3 * 10**5

    if not nx.has_path(contracted_adjacency_graph, from_end, to_start):
        return False

    return nx.dijkstra_path_length(contracted_adjacency_graph, from_end, to_start) < 3 * 10**5


def build_path_components(breakpoint_graph, max_matching):
    path_components = nx.Graph()

    for (node, data) in breakpoint_graph.nodes(data=True):
        path_components.add_node(node, **data)

    for (node_from, node_to) in max_matching:
        path_components.add_edge(node_from, node_to)

    for (node_from, data_from) in breakpoint_graph.nodes(data=True):
        for (node_to, data_to) in breakpoint_graph.nodes(data=True):
            if node_from == node_to:
                continue

            label_from, label_to = data_from["label"], data_to["label"]
            if label_from[:-1] == label_to[:-1]:
                path_components.add_edge(node_from, node_to)

    return path_components


def get_unused_edges(breakpoint_graph, max_matching):
    unused_edges = []

    for (node_from, node_to) in breakpoint_graph.edges():
        if (node_from, node_to) in max_matching:
            continue

        if (node_to, node_from) in max_matching:
            continue

        unused_edges.append((node_from, node_to))

    return unused_edges


def unite_cycles(graph, unused_edges):
    def is_cycle(start_node):
        for node in nx.dfs_preorder_nodes(graph, start_node):
            if graph.degree[node] != 2:
                return False

        return True

    def unite_components(node_from, node_to):
        node_from_label = node_labels[node_from]
        node_to_label = node_labels[node_to]

        for node in graph[node_from]:
            node_label = node_labels[node]
            if node_label[:-1] != node_from_label[:-1]:
                graph.remove_edge(node, node_from)
                break

        for node in graph[node_to]:
            node_label = node_labels[node]
            if node_label[:-1] != node_to_label[:-1]:
                graph.remove_edge(node, node_to)
                break

        graph.add_edge(node_from, node_to)

    number_united_components = 0
    node_labels = nx.get_node_attributes(graph, "label")
    for (node_from, node_to) in unused_edges:
        if nx.has_path(graph, node_from, node_to):
            continue

        component_from_is_cycle = is_cycle(node_from)
        component_to_is_cycle = is_cycle(node_to)

        if component_from_is_cycle and component_to_is_cycle:
            unite_components(node_from, node_to)
            number_united_components += 1
            continue

        if component_from_is_cycle and graph.degree[node_to] == 1:
            unite_components(node_from, node_to)
            number_united_components += 1
            continue

        if graph.degree[node_from] == 1 and component_to_is_cycle:
            unite_components(node_from, node_to)
            number_united_components += 1

    return number_united_components
