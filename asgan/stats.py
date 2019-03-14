import networkx as nx


def calc_stats(assembly_graph_query, alignment_blocks_query, full_paths_query,
               assembly_graph_target, alignment_blocks_target, full_paths_target, outdir):
    number_wcc_query = nx.number_weakly_connected_components(assembly_graph_query)
    number_wcc_target = nx.number_weakly_connected_components(assembly_graph_target)

    number_contigs_query = len(assembly_graph_query.edges())
    number_contigs_target = len(assembly_graph_target.edges())

    contig_lengths_query = [edge[2]["length"] for edge in assembly_graph_query.edges(data=True)]
    contig_lengths_target = [edge[2]["length"] for edge in assembly_graph_target.edges(data=True)]

    contigs_n50_query, _ = calc_nx(contig_lengths_query)
    contigs_n50_target, _ = calc_nx(contig_lengths_target)

    alignment_block_lengths_query = calc_alignment_block_lengths(alignment_blocks_query)
    alignment_block_lengths_target = calc_alignment_block_lengths(alignment_blocks_target)

    alignment_blocks_n50_query, _ = calc_nx(alignment_block_lengths_query)
    alignment_blocks_n50_target, _ = calc_nx(alignment_block_lengths_target)

    path_lengths_query = calc_path_lengths(full_paths_query)
    path_lengths_target = calc_path_lengths(full_paths_target)

    paths_n50_query, _ = calc_nx(path_lengths_query)
    paths_n50_target, _ = calc_nx(path_lengths_target)

    stats = dict()

    stats["number_wcc_query"] = number_wcc_query
    stats["number_wcc_target"] = number_wcc_target
    stats["number_contigs_query"] = number_contigs_query
    stats["number_contigs_target"] = number_contigs_target
    stats["contigs_n50_query"] = contigs_n50_query
    stats["contigs_n50_target"] = contigs_n50_target
    stats["alignment_blocks_n50_query"] = alignment_blocks_n50_query
    stats["alignment_blocks_n50_target"] = alignment_blocks_n50_target
    stats["paths_n50_query"] = paths_n50_query
    stats["paths_n50_target"] = paths_n50_target
    stats["number_paths"] = 2 * len(path_lengths_query)

    return stats


def calc_alignment_block_lengths(alignment_blocks):
    alignment_block_lengths = []

    for blocks in alignment_blocks.values():
        alignment_block_lengths.extend([block.length() for block in blocks])

    return alignment_block_lengths


def calc_path_lengths(paths):
    path_lengths = []

    for path in paths:
        path_lengths.append(calc_path_length(path[0]))

    return path_lengths


def calc_path_length(path):
    path_length = 0

    for i in range(len(path)):
        if i % 2 == 0:
            path_length += path[i].length()
        else:
            for (_, sequence_length) in path[i]:
                path_length += sequence_length

    if len(path) > 1 and path[0].signed_id() == path[-1].signed_id():
        path_length -= path[0].length()

    return path_length


def calc_nx(lengths, rate=0.5):
    total_length = sum(lengths)
    nx, lx, sum_length = 0, 0, 0

    for length in sorted(lengths, reverse=True):
        sum_length += length
        lx += 1
        if sum_length > rate * total_length:
            nx = length
            break

    return nx, lx
