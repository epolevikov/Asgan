import networkx as nx


def build_breakpoint_graph(graph_query, graph_target, num_aln_blocks):
    breakpoint_graph = nx.Graph()
    label2node = dict()

    for i in range(1, num_aln_blocks // 2 + 1):
        tail, head = "{}t".format(i), "{}h".format(i)

        breakpoint_graph.add_node(i, label=tail)
        breakpoint_graph.add_node(num_aln_blocks // 2 + i, label=head)

        label2node[tail] = i
        label2node[head] = num_aln_blocks // 2 + i

    block2edge_query = _build_block2edge_dict(graph_query)
    block2edge_target = _build_block2edge_dict(graph_target)
    contracted_graph_query = _build_contracted_graph(graph_query)
    contracted_graph_target = _build_contracted_graph(graph_target)
    signs = [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]

    for i in range(1, num_aln_blocks // 2 + 1):
        for j in range(i, num_aln_blocks // 2 + 1):
            for k, (sign1, sign2) in enumerate(signs):
                if i == j and k > 0:
                    break

                block_from_fwd, block_to_fwd = _signed_blocks(i, j,
                                                              sign1,
                                                              sign2)
                block_from_inv, block_to_inv = _signed_blocks(j, i,
                                                              _inv_sign(sign2),
                                                              _inv_sign(sign1))

                if not _check_block_pair(contracted_graph_query,
                                         block2edge_query,
                                         block_from_fwd, block_to_fwd):
                    continue

                if not _check_block_pair(contracted_graph_query,
                                         block2edge_query,
                                         block_from_inv, block_to_inv):
                    continue

                if not _check_block_pair(contracted_graph_target,
                                         block2edge_target,
                                         block_from_fwd, block_to_fwd):
                    continue

                if not _check_block_pair(contracted_graph_target,
                                         block2edge_target,
                                         block_from_inv, block_to_inv):
                    continue

                label_from = "{}{}".format(i, ["t", "h"][sign1 == "+"])
                label_to = "{}{}".format(j, ["h", "t"][sign2 == "+"])
                node_from = label2node[label_from]
                node_to = label2node[label_to]
                breakpoint_graph.add_edge(node_from, node_to)

    return breakpoint_graph


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


def _build_block2edge_dict(graph):
    block2edge = dict()
    for (node_from, node_to, data) in graph.edges(data=True):
        edge_name = data["name"]

        if edge_name.startswith("+") or edge_name.startswith("-"):
            block2edge[edge_name] = (node_from, node_to)

    return block2edge


def _build_contracted_graph(graph):
    contracted_graph = nx.MultiDiGraph()
    contracted_graph.add_nodes_from(graph.nodes)

    for (node_from, node_to, data) in graph.edges(data=True):
        if data["name"].startswith("+") or data["name"].startswith("-"):
            continue

        contracted_graph.add_edge(node_from, node_to, **data)

    return contracted_graph


def _signed_blocks(block_from, block_to, sign_from, sign_to):
    block_from_signed = "{}{}".format(sign_from, block_from)
    block_to_signed = "{}{}".format(sign_to, block_to)
    return block_from_signed, block_to_signed


def _check_block_pair(contracted_graph, block2edge, block_from, block_to):
    (node_from1, node_to1) = block2edge[block_from]
    (node_from2, node_to2) = block2edge[block_to]
    return nx.has_path(contracted_graph, node_to1, node_from2)


def _inv_sign(sign):
    return ["+", "-"][sign == "+"]
