

def output_blocks_info(out_dir, aln_blocks_query, aln_blocks_target):
    def _output_blocks_info(file, aln_blocks):
        for seq, blocks in aln_blocks.items():
            blocks.sort(key=lambda block: block.start)
            for block in blocks:
                f.write(block.signed_id() + " ")
            f.write("| seq={} len={}\n".format(
                    blocks[0].seq_name,
                    pretty_number(blocks[0].seq_length)))

    with open("{}/alignment_blocks.txt".format(out_dir), "w") as f:
        f.write("Query blocks:\n")
        _output_blocks_info(f, aln_blocks_query)
        f.write("\nTarget blocks:\n")
        _output_blocks_info(f, aln_blocks_target)


def alignment_graph_save_dot(graph, outdir, outfile):
    def get_edge_color(edge_name):
        if edge_name.startswith("+") or edge_name.startswith("-"):
            return "green"

        return "black"

    with open("{}/{}.gv".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point]\n")
        f.write("  edge [penwidth=5]\n")

        for (node_from, node_to, data) in graph.edges(data=True):
            f.write("  {} -> {} [label=\"{}\", color={}]\n".format(
                node_from, node_to, data["name"],
                get_edge_color(data["name"])))

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
