import networkx as nx


def build_breakpoint_graph(alignment_graph_query, alignment_blocks_query,
                           alignment_graph_target, alignment_blocks_target):
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

    id2block_query = build_id2block_dict(alignment_blocks_query)
    id2block_target = build_id2block_dict(alignment_blocks_target)

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

                if not _check_block_pair(contracted_alignment_graph_query,
                                         block2edge_query, id2block_query,
                                         block_from_fwd, block_to_fwd):
                    continue

                if not _check_block_pair(contracted_alignment_graph_query,
                                         block2edge_query, id2block_query,
                                         block_from_inv, block_to_inv):
                    continue

                if not _check_block_pair(contracted_alignment_graph_target,
                                         block2edge_target, id2block_target,
                                         block_from_fwd, block_to_fwd):
                    continue

                if not _check_block_pair(contracted_alignment_graph_target,
                                         block2edge_target, id2block_target,
                                         block_from_inv, block_to_inv):
                    continue

                label_from = "{}{}".format(i, ["t", "h"][sign1 == "+"])
                label_to = "{}{}".format(j, ["h", "t"][sign2 == "+"])
                node_from = label2node[label_from]
                node_to = label2node[label_to]
                breakpoint_graph.add_edge(node_from, node_to)

    return breakpoint_graph


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


def unite_cycles(graph, unused_edges, args):
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

    # with open("{}/united_components.txt".format(args.out_dir), "w") as f:
    #    f.write("number united components: {}\n".format(number_united_components))


def count_number_alignment_blocks(alignment_graph):
    number_alignment_blocks = 0

    for (node_from, node_to, data) in alignment_graph.edges(data=True):
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


def build_id2block_dict(alignment_blocks):
    id2block = dict()

    for blocks in alignment_blocks.values():
        for block in blocks:
            id2block[block.signed_id()] = block

    return id2block


def build_contracted_alignment_graph(alignment_graph):
    contracted_alignment_graph = nx.MultiDiGraph()

    for node, data in alignment_graph.nodes(data=True):
        contracted_alignment_graph.add_node(node, **data)

    for node_from, node_to, data in alignment_graph.edges(data=True):
        if data["name"].startswith("+") or data["name"].startswith("-"):
            continue

        contracted_alignment_graph.add_edge(node_from, node_to,
                                            name=data["name"],
                                            weight=data["length"])

    return contracted_alignment_graph


def _make_signed_blocks(block_id_from, block_id_to, sign_from, sign_to):
    block_from_signed = "{}{}".format(sign_from, block_id_from)
    block_to_signed = "{}{}".format(sign_to, block_id_to)
    return block_from_signed, block_to_signed


def _check_block_pair(contracted_alignment_graph, block2edge, id2block, block_id_from, block_id_to):
    (from_start, from_end) = block2edge[block_id_from]
    (to_start, to_end) = block2edge[block_id_to]

    if from_end == to_start:
        dist = contracted_alignment_graph.nodes[from_end].get("distance")

        if dist is None:
            block_from = id2block[block_id_from]
            block_to = id2block[block_id_to]

            dist = block_to.start + (block_from.seq_length - block_from.end)

        return dist < 3 * 10**5

    if not nx.has_path(contracted_alignment_graph, from_end, to_start):
        return False

    return nx.dijkstra_path_length(contracted_alignment_graph, from_end, to_start) < 3 * 10**5


def _inv_sign(sign):
    return ["+", "-"][sign == "+"]
