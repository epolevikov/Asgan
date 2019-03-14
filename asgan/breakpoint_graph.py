import networkx as nx


def build_breakpoint_graph(alignment_graph_query, alignment_graph_target):
    number_alignment_blocks = count_number_alignment_blocks(alignment_graph_query)
    breakpoint_graph = nx.Graph()
    label2node = dict()

    for i in range(1, number_alignment_blocks // 2 + 1):
        tail, head = "{}t".format(i), "{}h".format(i)

        breakpoint_graph.add_node(i, label=tail)
        breakpoint_graph.add_node(number_alignment_blocks // 2 + i, label=head)

        label2node[tail] = i
        label2node[head] = number_alignment_blocks // 2 + i

    block2edge_query = build_block2edge_dict(alignment_graph_query)
    block2edge_target = build_block2edge_dict(alignment_graph_target)

    contracted_alignment_graph_query = build_contracted_alignment_graph(alignment_graph_query)
    contracted_alignment_graph_target = build_contracted_alignment_graph(alignment_graph_target)

    signs = [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]

    for i in range(1, number_alignment_blocks // 2 + 1):
        for j in range(i, number_alignment_blocks // 2 + 1):
            for k, (sign1, sign2) in enumerate(signs):
                if i == j and k > 0:
                    break

                block_from_fwd, block_to_fwd = _make_signed_blocks(i, j, sign1, sign2)
                block_from_inv, block_to_inv = _make_signed_blocks(j, i, _inv_sign(sign2), _inv_sign(sign1))

                if not _check_block_pair(contracted_alignment_graph_query, block2edge_query,
                                         block_from_fwd, block_to_fwd):
                    continue

                if not _check_block_pair(contracted_alignment_graph_query, block2edge_query,
                                         block_from_inv, block_to_inv):
                    continue

                if not _check_block_pair(contracted_alignment_graph_target, block2edge_target,
                                         block_from_fwd, block_to_fwd):
                    continue

                if not _check_block_pair(contracted_alignment_graph_target, block2edge_target,
                                         block_from_inv, block_to_inv):
                    continue

                label_from = "{}{}".format(i, ["t", "h"][sign1 == "+"])
                label_to = "{}{}".format(j, ["h", "t"][sign2 == "+"])
                node_from = label2node[label_from]
                node_to = label2node[label_to]
                breakpoint_graph.add_edge(node_from, node_to)

    return breakpoint_graph


def count_number_alignment_blocks(alignmnet_graph):
    number_alignment_blocks = 0

    for (node_from, node_to, data) in alignmnet_graph.edges(data=True):
        if data["name"].startswith("+") or data["name"].startswith("-"):
            number_alignment_blocks += 1

    return number_alignment_blocks


def get_unused_edges(breakpoint_graph, max_matching):
    unused_edges = []

    for (node_from, node_to) in breakpoint_graph.edges():
        if (node_from, node_to) in max_matching:
            continue

        if (node_to, node_from) in max_matching:
            continue

        unused_edges.append((node_from, node_to))

    return unused_edges


def build_block2edge_dict(alignment_graph):
    block2edge = dict()

    for (node_from, node_to, data) in alignment_graph.edges(data=True):
        edge_name = data["name"]

        if edge_name.startswith("+") or edge_name.startswith("-"):
            block2edge[edge_name] = (node_from, node_to)

    return block2edge


def build_contracted_alignment_graph(alignment_graph):
    contracted_alignment_graph = nx.MultiDiGraph()

    for node, data in alignment_graph.nodes(data=True):
        contracted_alignment_graph.add_node(node, **data)

    for node_from, node_to, data in alignment_graph.edges(data=True):
        if data["name"].startswith("+") or data["name"].startswith("-"):
            continue

        contracted_alignment_graph.add_edge(node_from, node_to, name=data["name"], weight=data["length"])

    return contracted_alignment_graph


def _make_signed_blocks(block_id_from, block_id_to, sign_from, sign_to):
    block_from_signed = "{}{}".format(sign_from, block_id_from)
    block_to_signed = "{}{}".format(sign_to, block_id_to)
    return block_from_signed, block_to_signed


def _check_block_pair(contracted_alignment_graph, block2edge, block_id_from, block_id_to):
    (from_start, from_end) = block2edge[block_id_from]
    (to_start, to_end) = block2edge[block_id_to]
    return nx.has_path(contracted_alignment_graph, from_end, to_start)


def _inv_sign(sign):
    return ["+", "-"][sign == "+"]
