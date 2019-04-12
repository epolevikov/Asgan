

def assembly_graph_save_dot(graph, out_dir, out_file):
    with open("{}/{}".format(out_dir, out_file), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point]\n")
        f.write("  edge [penwidth=5, color=green, fontsize=20]\n")

        for (node_from, node_to, data) in graph.edges(data=True):
            edge_color = ["green", "black"][data["is_repeat"]]
            f.write("  {} -> {} [label=\"{}\", color=\"{}\"]\n".format(
                node_from, node_to, data["name"], edge_color))

        f.write("}\n")


def save_raw_hits(raw_hits, out_dir, out_file):
    raw_hits.sort(key=lambda hit: (hit.query_name, hit.query_start))
    out_file = "{}/{}".format(out_dir, out_file)

    with open(out_file, "w") as f:
        for raw_hit in raw_hits:
            f.write(str(raw_hit) + "\n")


def output_blocks_info(synteny_blocks_query, synteny_blocks_target, out_dir, out_file):
    def _output_blocks_info(synteny_blocks):
        for seq, blocks in synteny_blocks.items():
            blocks.sort(key=lambda block: block.start)

            sequence_name = blocks[0].sequence_name
            sequence_length = pretty_number(blocks[0].sequence_length)
            f.write("seq={} len={}\n".format(sequence_name, sequence_length))

            for block in blocks:
                f.write(str(block) + "\n")

            f.write("\n")

    with open("{}/{}".format(out_dir, out_file), "w") as f:
        f.write("Query blocks:\n")
        _output_blocks_info(synteny_blocks_query)

        f.write("\nTarget blocks:\n")
        _output_blocks_info(synteny_blocks_target)


def adjacency_graph_save_dot(graph, out_dir, out_file, block_attributes=None):
    def get_edge_color(edge_name):
        if block_attributes is None:
            if edge_name.startswith("+") or edge_name.startswith("-"):
                return "green"
            else:
                return "black"
        else:
            edge_color = block_attributes["color"].get(edge_name)

        if edge_color is None:
            return "black"

        return edge_color

    def get_edge_label(edge_name):
        if edge_name.startswith("+") or edge_name.startswith("-"):
            return edge_name

        return ""

    def get_edge_style(edge_name):
        if block_attributes is None:
            return "solid"
        else:
            edge_style = block_attributes["style"].get(edge_name)

            if edge_style is None:
                return "solid"

            return edge_style

    with open("{}/{}".format(out_dir, out_file), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point, width=0.06]\n")
        f.write("  edge [fontsize=20]\n")
        f.write("  graph[center=true, margin=0.5, ")
        f.write("nodesep=0.45, ranksep=0.35]\n")

        # for node, data in graph.nodes(data=True):
        #    dist = data.get("distance")
        #    f.write("  {} [label=\"{}\"]\n".format(node, dist))

        for node_from, node_to, data in graph.edges(data=True):
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


def breakpoint_graph_save_dot(graph, max_matching, out_dir, out_file):
    def get_edge_color(node_from, node_to):
        if (node_from, node_to) in max_matching:
            return "green"

        if (node_to, node_from) in max_matching:
            return "green"

        return "black"

    with open("{}/{}".format(out_dir, out_file), "w") as f:
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


def path_components_save_dot(paths, unused_edges, out_dir):
    def save(out_file, with_unused_edges=False):
        with open("{}/{}.gv".format(out_dir, out_file), "w") as f:
            f.write("graph {\n")
            f.write("  node [fontsize=20]\n")
            f.write("  edge [penwidth=5]\n")

            for (node, data) in paths.nodes(data=True):
                f.write("  {} [label=\"{}\"]\n".format(node, data["label"]))

            for (node_from, node_to) in paths.edges():
                f.write("  {} -- {} [color=green] \n".format(node_from, node_to))

            if with_unused_edges:
                for (node_from, node_to) in unused_edges:
                    f.write("  {} -- {} [style=dashed]\n".format(node_from, node_to))

            f.write("}\n")

    save("paths")
    save("paths_with_unused_edges", with_unused_edges=True)


def save_path_sequences(paths_query, paths_target, out_dir):
    def write_path(path, f):
        for i in range(len(path)):
            if i % 2 == 0:
                block_id = path[i].signed_id()
                block_length = path[i].length()
                block_length = pretty_number(block_length)
                f.write("{} [{}]\n".format(block_id, block_length))
            else:
                for (sequence_name, sequence_length) in path[i]:
                    sequence_length = pretty_number(sequence_length)
                    f.write("  {} [{}]\n".format(sequence_name, sequence_length))

    with open("{}/paths.txt".format(out_dir), "w") as f:
        for i in range(len(paths_query)):
            f.write("query path:\n")
            write_path(paths_query[i], f)
            f.write("\ntarget path:\n")
            write_path(paths_target[i], f)
            f.write("------------\n\n")


def output_stats(stats, out_dir):
    with open("{}/stats.txt".format(out_dir), "w") as f:
        # header
        f.write("\tQuery       \tTarget\n")

        # wcc
        number_min_width = 12
        number_wcc_query = pretty_number(stats["number_wcc_query"], number_min_width)
        number_wcc_target = pretty_number(stats["number_wcc_target"], number_min_width)
        f.write("cc\t{}\t{}\n\n".format(number_wcc_query, number_wcc_target))

        # sequences
        number_sequences_query = pretty_number(stats["number_sequences_query"], number_min_width)
        number_sequences_target = pretty_number(stats["number_sequences_target"], number_min_width)
        sequences_total_length_query = pretty_number(stats["sequences_total_length_query"], number_min_width)
        sequences_total_length_target = pretty_number(stats["sequences_total_length_target"], number_min_width)
        sequences_n50_query = pretty_number(stats["sequences_n50_query"], number_min_width)
        sequences_n50_target = pretty_number(stats["sequences_n50_target"], number_min_width)
        sequences_l50_query = pretty_number(stats["sequences_l50_query"], number_min_width)
        sequences_l50_target = pretty_number(stats["sequences_l50_target"], number_min_width)

        f.write("seqs\t{}\t{}\n".format(number_sequences_query, number_sequences_target))
        f.write("tlen\t{}\t{}\n".format(sequences_total_length_query, sequences_total_length_target))
        f.write("N50\t{}\t{}\n".format(sequences_n50_query, sequences_n50_target))
        f.write("L50\t{}\t{}\n\n".format(sequences_l50_query, sequences_l50_target))

        # blocks
        number_blocks_query = pretty_number(stats["number_blocks"], number_min_width)
        number_blocks_target = pretty_number(stats["number_blocks"], number_min_width)
        blocks_total_length_query = pretty_number(stats["blocks_total_length_query"], number_min_width)
        blocks_total_length_target = pretty_number(stats["blocks_total_length_target"], number_min_width)
        blocks_n50_query = pretty_number(stats["blocks_n50_query"], number_min_width)
        blocks_n50_target = pretty_number(stats["blocks_n50_target"], number_min_width)
        blocks_l50_query = pretty_number(stats["blocks_l50_query"], number_min_width)
        blocks_l50_target = pretty_number(stats["blocks_l50_target"], number_min_width)

        f.write("blocks\t{}\t{}\n".format(number_blocks_query, number_blocks_target))
        f.write("tlen\t{}\t{}\n".format(blocks_total_length_query, blocks_total_length_target))
        f.write("N50\t{}\t{}\n".format(blocks_n50_query, blocks_n50_target))
        f.write("L50\t{}\t{}\n\n".format(blocks_l50_query, blocks_l50_target))

        # paths
        number_paths_query = pretty_number(stats["number_paths"], number_min_width)
        number_paths_target = pretty_number(stats["number_paths"], number_min_width)
        paths_total_length_query = pretty_number(stats["paths_total_length_query"], number_min_width)
        paths_total_length_target = pretty_number(stats["paths_total_length_target"], number_min_width)
        paths_n50_query = pretty_number(stats["paths_n50_query"], number_min_width)
        paths_n50_target = pretty_number(stats["paths_n50_target"], number_min_width)
        paths_l50_query = pretty_number(stats["paths_l50_query"], number_min_width)
        paths_l50_target = pretty_number(stats["paths_l50_target"], number_min_width)

        f.write("paths\t{}\t{}\n".format(number_paths_query, number_paths_target))
        f.write("tlen\t{}\t{}\n".format(paths_total_length_query, paths_total_length_target))
        f.write("N50\t{}\t{}\n".format(paths_n50_query, paths_n50_target))
        f.write("L50\t{}\t{}\n\n".format(paths_l50_query, paths_l50_target))

        # link types
        f.write("link types: {} {} {} {}".format(*stats["link_types"]))


def pretty_number(number, min_width=None, fill=" "):
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

    if min_width is not None and len(number) < min_width:
        number += " " * (min_width - len(number))

    return number


'''
from asgan.stats import calc_path_length


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

    with open("{}/{}".format(outdir, outfile), "w") as f:
        f.write("digraph {\n")
        f.write("  node [shape=point, width=0.06]\n")
        f.write("  edge [fontsize=20]\n")
        f.write("  graph[center=true, margin=0.5, ")
        f.write("nodesep=0.45, ranksep=0.35]\n")

        # for node, data in graph.nodes(data=True):
        #    dist = data.get("distance")
        #    f.write("  {} [label=\"{}\"]\n".format(node, dist))

        for node_from, node_to, data in graph.edges(data=True):
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


def save_full_paths(paths_query, paths_target, outdir):
    def write_path(path, f):
        path_length = calc_path_length(path)

        for i in range(len(path)):
            if i % 2 == 0:
                f.write("{} [{}]\n".format(path[i].signed_id(), pretty_number(path[i].length())))
            else:
                for (sequence_name, sequence_length) in path[i]:
                    f.write("  {} [{}]\n".format(sequence_name, pretty_number(sequence_length)))

        f.write("length={}\n\n".format(pretty_number(path_length)))

    with open("{}/paths.txt".format(outdir), "w") as f:
        for i in range(len(paths_query)):
            f.write("query forward:\n")
            write_path(paths_query[i][0], f)
            f.write("target forward:\n")
            write_path(paths_target[i][0], f)

            f.write("query reverse:\n")
            write_path(paths_query[i][1], f)
            f.write("target reverse:\n")
            write_path(paths_target[i][1], f)
            f.write("------------\n\n")


def output_stats(stats, outdir):
    min_length = 12
    with open("{}/stats.txt".format(outdir), "w") as f:
        f.write("\tQuery       \tTarget\n")
        f.write("CC\t{}\t{}\n\n".format(pretty_number(stats["number_wcc_query"], min_length=min_length),
                                        pretty_number(stats["number_wcc_target"], min_length=min_length)))

        f.write("contigs\t{}\t{}\n".format(pretty_number(stats["number_contigs_query"], min_length=min_length),
                                           pretty_number(stats["number_contigs_target"], min_length=min_length)))
        f.write("tlen\t{}\t{}\n".format(pretty_number(stats["contigs_total_length_query"], min_length=min_length),
                                        pretty_number(stats["contigs_total_length_target"], min_length=min_length)))
        f.write("N50\t{}\t{}\n".format(pretty_number(stats["contigs_n50_query"], min_length=min_length),
                                       pretty_number(stats["contigs_n50_target"], min_length=min_length)))
        f.write("L50\t{}\t{}\n\n".format(pretty_number(stats["contigs_l50_query"], min_length=min_length),
                                         pretty_number(stats["contigs_l50_target"], min_length=min_length)))

        f.write("blocks\t{}\t{}\n".format(pretty_number(stats["number_alignment_blocks"], min_length=min_length),
                                          pretty_number(stats["number_alignment_blocks"], min_length=min_length)))
        f.write("tlen\t{}\t{}\n".format(pretty_number(stats["alignment_blocks_total_length_query"],
                                                      min_length=min_length),
                                        pretty_number(stats["alignment_blocks_total_length_target"],
                                                      min_length=min_length)))
        f.write("N50\t{}\t{}\n".format(pretty_number(stats["alignment_blocks_n50_query"], min_length=min_length),
                                       pretty_number(stats["alignment_blocks_n50_target"], min_length=min_length)))
        f.write("L50\t{}\t{}\n\n".format(pretty_number(stats["alignment_blocks_l50_query"], min_length=min_length),
                                         pretty_number(stats["alignment_blocks_l50_target"], min_length=min_length)))

        f.write("s-paths\t{}\t{}\n".format(pretty_number(stats["number_single_paths_query"], min_length=min_length),
                                           pretty_number(stats["number_single_paths_target"], min_length=min_length)))
        f.write("tlen\t{}\t{}\n".format(pretty_number(stats["single_paths_total_length_query"], min_length=min_length),
                                        pretty_number(stats["single_paths_total_length_target"],
                                        min_length=min_length)))
        f.write("N50\t{}\t{}\n".format(pretty_number(stats["single_paths_n50_query"], min_length=min_length),
                                       pretty_number(stats["single_paths_n50_target"], min_length=min_length)))
        f.write("L50\t{}\t{}\n\n".format(pretty_number(stats["single_paths_l50_query"], min_length=min_length),
                                         pretty_number(stats["single_paths_l50_target"], min_length=min_length)))

        f.write("u-paths\t{}\t{}\n".format(pretty_number(stats["number_paths"], min_length=min_length),
                                           pretty_number(stats["number_paths"], min_length=min_length)))
        f.write("tlen\t{}\t{}\n".format(pretty_number(stats["paths_total_length_query"], min_length=min_length),
                                        pretty_number(stats["paths_total_length_target"], min_length=min_length)))
        f.write("N50\t{}\t{}\n".format(pretty_number(stats["paths_n50_query"], min_length=min_length),
                                       pretty_number(stats["paths_n50_target"], min_length=min_length)))
        f.write("L50\t{}\t{}\n".format(pretty_number(stats["paths_l50_query"], min_length=min_length),
                                       pretty_number(stats["paths_l50_target"], min_length=min_length)))


        f.write("\nlink types: {} {} {} {}".format(*stats["link_types"]))


def pretty_number(number, min_width=None, fill=" "):
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

    if min_length is not None and len(number) < min_length:
        number += " " * (min_length - len(number))

    return number
'''
