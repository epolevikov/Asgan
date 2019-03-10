

def output_blocks_info(aln_blocks_query, aln_blocks_target, out_dir):
    def _output_blocks_info(file, aln_blocks):
        for seq, blocks in aln_blocks.items():
            blocks.sort(key=lambda block: block.start)

            seq_name = blocks[0].seq_name
            seq_len = pretty_number(blocks[0].seq_length)
            f.write("seq={} len={}\n".format(seq_name, seq_len))

            for block in blocks:
                f.write(str(block) + "\n")

            f.write("\n")

    with open("{}/alignment_blocks.txt".format(out_dir), "w") as f:
        f.write("Query blocks:\n")
        _output_blocks_info(f, aln_blocks_query)

        f.write("\nTarget blocks:\n")
        _output_blocks_info(f, aln_blocks_target)


def assembly_graph_save_dot(graph, outfile, outdir):
    with open("{}/{}.gv".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point]\n")
        f.write("  edge [penwidth=5, color=green, fontsize=20]\n")

        for (node_from, node_to, data) in graph.edges(data=True):
            f.write("  {} -> {} [label=\"{}\"]\n".format(node_from, node_to, data["name"]))

        f.write("}\n")


def alignment_graph_save_dot(graph, outfile, block_colors, block_styles, outdir):
    def get_edge_color(edge_name):
        edge_color = block_colors.get(edge_name)

        if edge_color is not None:
            return edge_color

        return "black"

    def get_edge_label(edge_name):
        if edge_name.startswith("+") or edge_name.startswith("-"):
            return edge_name

        return ""

    def get_edge_style(edge_name):
        edge_style = block_styles.get(edge_name)

        if edge_style is not None:
            return edge_style

        return "solid"

    with open("{}/{}.gv".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point, width=0.06]\n")
        f.write("  edge [fontsize=20]\n")
        f.write("  graph[center=true, margin=0.5, ")
        f.write("nodesep=0.45, ranksep=0.35]\n")

        for (node_from, node_to, data) in graph.edges(data=True):
            edge_color = get_edge_color(data["name"])
            edge_label = get_edge_label(data["name"])
            edge_style = get_edge_style(data["name"])
            edge_penwidth = 3 if edge_color == "black" else 6

            f.write("  {} -> {} [".format(node_from, node_to))
            f.write("label=\"{}\", ".format(edge_label))
            f.write("color=\"{}\", ".format(edge_color))
            f.write("style=\"{}\", ".format(edge_style))
            f.write("penwidth=\"{}\"]\n".format(edge_penwidth))

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
        f.write("  node [fontsize=20]\n")

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
            f.write("  node [fontsize=20]\n")
            f.write("  edge [penwidth=5]\n")

            for node, data in paths.nodes(data=True):
                f.write("  {} [label=\"{}\"]\n".format(node, data["label"]))

            for (node_from, node_to) in paths.edges():
                f.write("  {} -- {} [color=green] \n".format(node_from, node_to))

            if with_unused_edges:
                for (node_from, node_to) in unused_edges:
                    f.write("  {} -- {} [style=dashed]\n".format(node_from, node_to))

            f.write("}\n")

    save("paths")
    save("paths_with_unused_edges", with_unused_edges=True)


def output_stats(stats, outdir):
    with open("{}/stats.txt".format(outdir), "w") as f:
        f.write("\twcc\tcontigs\tpaths\n")
        f.write("Query\t{}\t{}\t{}\n".format(stats["query_num_wcc"],
                                             stats["query_num_contigs"],
                                             2 * stats["num_paths"]))
        f.write("Target\t{}\t{}\t{}\n".format(stats["target_num_wcc"],
                                              stats["target_num_contigs"],
                                              2 * stats["num_paths"]))


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
