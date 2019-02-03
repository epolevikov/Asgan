

def output_blocks_info(out_dir, aln_blocks_query, aln_blocks_target):
    def _output_blocks_info(file, aln_blocks):
        for seq, blocks in aln_blocks.items():
            blocks.sort(key=lambda block: block.start)
            '''
            for block in blocks:
                f.write(block.signed_id() + " ")
            '''
            f.write("seq={} len={}\n".format(
                    blocks[0].seq_name,
                    pretty_number(blocks[0].seq_length)))
            for block in blocks:
                f.write(str(block) + "\n")
            f.write("\n")

    with open("{}/alignment_blocks.txt".format(out_dir), "w") as f:
        f.write("Query blocks:\n")
        _output_blocks_info(f, aln_blocks_query)
        f.write("\nTarget blocks:\n")
        _output_blocks_info(f, aln_blocks_target)


def assembly_graph_save_dot(graph, outdir, outfile):
    with open("{}/{}.gv".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point]\n")
        f.write("  edge [penwidth=5]\n")

        for (node_from, node_to, data) in graph.edges(data=True):
            f.write("  {} -> {} [label=\"{}\"]\n".format(
                node_from, node_to, data["name"]))

        f.write("}\n")


'''
def alignment_graph_save_dot(graph, label, prefix, block_colors, file):
    def get_edge_color(edge_name):
        edge_color = block_colors.get(edge_name)
        return [edge_color, "black"][edge_color is None]

    def get_edge_label(edge_name):
        if edge_name.startswith("+") or edge_name.startswith("-"):
            return edge_name

        return ""

    file.write("  subgraph cluster_{} ".format(prefix) + " {\n")

    for (node_from, node_to, data) in graph.edges(data=True):
        edge_color = get_edge_color(data["name"])
        edge_label = get_edge_label(data["name"])
        file.write("    {0}{1} -> {0}{2} [label=\"{3}\", color={4}]\n".format(
            prefix, node_from, node_to, edge_label, edge_color))

    file.write("    label={}\n".format(label))
    file.write("  }\n")


def alignment_graphs_save_dot(graph_query, graph_target, block_colors, outdir,
                              outfile):
    with open("{}/{}".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point]\n")
        f.write("  edge [penwidth=5]\n")
        f.write("  graph[center=true, margin=0.2,")
        f.write("nodesep=0.35, ranksep=0.35]\n")

        alignment_graph_save_dot(graph_query, "Query", "qry",
                                 block_colors, f)
        alignment_graph_save_dot(graph_target, "Target", "trg",
                                 block_colors, f)

        f.write("}\n")
'''


def alignment_graph_save_dot(graph, block_colors, outdir, outfile):
    def get_edge_color(edge_name):
        edge_color = block_colors.get(edge_name)
        return [edge_color, "black"][edge_color is None]

    def get_edge_label(edge_name):
        if edge_name.startswith("+") or edge_name.startswith("-"):
            return edge_name

        return ""

    with open("{}/{}.gv".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point]\n")
        f.write("  edge [penwidth=5]\n")
        f.write("  graph[center=true, margin=0.2,")
        f.write("nodesep=0.35, ranksep=0.35]\n")

        for (node_from, node_to, data) in graph.edges(data=True):
            edge_color = get_edge_color(data["name"])
            edge_label = get_edge_label(data["name"])
            f.write("  {} -> {} [label=\"{}\", color={}]\n".format(
                node_from, node_to, edge_label, edge_color))

        f.write("}\n")


def breakpoint_graph_save_dot(graph, max_matching, outdir):
    def get_edge_color(node_from, node_to):
        if (node_from, node_to) in max_matching:
            return "green"

        if (node_to, node_from) in max_matching:
            return "green"

        return "black"

    with open("{}/{}.gv".format(outdir, "breakpoint_graph"), "w") as f:
        f.write("graph {\n")
        f.write("  edge [penwidth=5]\n")

        for node, data in graph.nodes(data=True):
            f.write("  {} [label=\"{}\"]\n".format(node, data["label"]))

        for (node_from, node_to) in graph.edges():
            color = get_edge_color(node_from, node_to)
            f.write("  {} -- {} [color={}]\n".format(
                node_from, node_to, color))

        f.write("}\n")


def paths_graph_save_dot(paths, unused_edges, out_dir):
    def save(filename, with_unused_edges=False):
        with open("{}/{}.gv".format(out_dir, filename), "w") as f:
            f.write("graph {\n")
            f.write("  edge [penwidth=5]\n")

            for node, data in paths.nodes(data=True):
                f.write("  {} [label=\"{}\"]\n".format(node, data["label"]))

            for (node_from, node_to) in paths.edges():
                f.write("  {} -- {} [color=green] \n".format(
                    node_from, node_to))

            if with_unused_edges:
                for (node_from, node_to) in unused_edges:
                    f.write("  {} -- {} [style=dashed]\n".format(
                        node_from, node_to))

            f.write("}\n")

    save("paths")
    save("paths_with_unused_edges", with_unused_edges=True)


def output_stats(stats, outdir):
    with open("{}/stats.txt".format(outdir), "w") as f:
        f.write("# query wcc: {}\n".format(stats.query_num_wcc))
        f.write("# query contigs: {}\n\n".format(stats.query_num_contigs))
        f.write("# target wcc: {}\n".format(stats.target_num_wcc))
        f.write("# target contigs: {}\n\n".format(stats.target_num_contigs))
        f.write("# paths: {}\n".format(2 * len(stats.paths)))


def pretty_number(number):
    digits = []
    number = str(number)
    is_neg = False

    if number.startswith("-"):
        is_neg = True
        number = number[1:]

    start = int(round(len(number) % 3, 0))

    if start != 0:
        digits.append(number[:start])
    digits.extend([number[i:i+3] for i in range(start, len(number), 3)])

    number = "\'".join(digits)
    if is_neg:
        number = "-" + number

    return number
