import networkx as nx


def calc_stats(assembly_graph_query, alignment_blocks_query, full_paths_query,
               assembly_graph_target, alignment_blocks_target, full_paths_target,
               alignment_block_paths, outdir):
    number_wcc_query = nx.number_weakly_connected_components(assembly_graph_query)
    number_wcc_target = nx.number_weakly_connected_components(assembly_graph_target)

    # contigs
    contig_lengths_query = get_contig_lengths(assembly_graph_query)
    contig_lengths_target = get_contig_lengths(assembly_graph_target)

    number_contigs_query = len(contig_lengths_query)
    number_contigs_target = len(contig_lengths_target)

    contigs_total_length_query = sum(contig_lengths_query)
    contigs_total_length_target = sum(contig_lengths_target)

    contigs_n50_query, contigs_l50_query = calc_nx(contig_lengths_query)
    contigs_n50_target, contigs_l50_target = calc_nx(contig_lengths_target)

    # alignment blocks
    alignment_block_lengths_query = get_alignment_block_lengths(alignment_blocks_query)
    alignment_block_lengths_target = get_alignment_block_lengths(alignment_blocks_target)

    alignment_blocks_total_length_query = sum(alignment_block_lengths_query)
    alignment_blocks_total_length_target = sum(alignment_block_lengths_target)

    number_alignment_blocks = len(alignment_block_lengths_query)

    alignment_blocks_n50_query, alignment_blocks_l50_query = calc_nx(alignment_block_lengths_query)
    alignment_blocks_n50_target, alignment_blocks_l50_target = calc_nx(alignment_block_lengths_target)

    # paths
    path_lengths_query = get_path_lengths(full_paths_query)
    path_lengths_target = get_path_lengths(full_paths_target)

    paths_total_length_query = sum(path_lengths_query)
    paths_total_length_target = sum(path_lengths_target)

    number_paths = len(path_lengths_query)

    paths_n50_query, paths_l50_query = calc_nx(path_lengths_query)
    paths_n50_target, paths_l50_target = calc_nx(path_lengths_target)

    # link types:

    link_types = calc_link_types(alignment_block_paths, alignment_blocks_query, alignment_blocks_target)

    stats = dict()

    stats["link_types"] = link_types

    stats["number_wcc_query"] = number_wcc_query
    stats["number_wcc_target"] = number_wcc_target

    stats["number_contigs_query"] = number_contigs_query
    stats["number_contigs_target"] = number_contigs_target
    stats["contigs_n50_query"] = contigs_n50_query
    stats["contigs_n50_target"] = contigs_n50_target
    stats["contigs_l50_query"] = contigs_l50_query
    stats["contigs_l50_target"] = contigs_l50_target
    stats["contigs_total_length_query"] = contigs_total_length_query
    stats["contigs_total_length_target"] = contigs_total_length_target

    stats["number_alignment_blocks"] = number_alignment_blocks
    stats["alignment_blocks_n50_query"] = alignment_blocks_n50_query
    stats["alignment_blocks_n50_target"] = alignment_blocks_n50_target
    stats["alignment_blocks_l50_query"] = alignment_blocks_l50_query
    stats["alignment_blocks_l50_target"] = alignment_blocks_l50_target
    stats["alignment_blocks_total_length_query"] = alignment_blocks_total_length_query
    stats["alignment_blocks_total_length_target"] = alignment_blocks_total_length_target

    stats["number_paths"] = number_paths
    stats["paths_n50_query"] = paths_n50_query
    stats["paths_n50_target"] = paths_n50_target
    stats["paths_l50_query"] = paths_l50_query
    stats["paths_l50_target"] = paths_l50_target
    stats["paths_total_length_query"] = paths_total_length_query
    stats["paths_total_length_target"] = paths_total_length_target

    return stats


def calc_link_types(alignment_block_paths, alignment_blocks_query, alignment_blocks_target):
    def get_link_type(block_from, block_to, id2block_query, id2block_target):
        link_types = [[0, 2], [1, 3]]

        seq_from_query = id2block_query[block_from].seq_name
        seq_to_query = id2block_query[block_to].seq_name

        seq_from_target = id2block_target[block_from].seq_name
        seq_to_target = id2block_target[block_to].seq_name

        row = seq_from_query != seq_to_query
        col = seq_from_target != seq_to_target

        # print(block_from, block_to)
        # print(seq_from_query, seq_to_query)
        # print(seq_from_target, seq_to_target)

        return link_types[row][col]

    id2block_query = build_id2block_dict(alignment_blocks_query)
    id2block_target = build_id2block_dict(alignment_blocks_target)
    link_types = [0, 0, 0, 0]

    for path in alignment_block_paths:
        if len(path[0]) > 1:
            for i in range(len(path[0]) - 1):
                block_from, block_to = path[0][i], path[0][i + 1]
                link_type = get_link_type(block_from, block_to, id2block_query, id2block_target)
                link_types[link_type] += 1

    return link_types


def build_id2block_dict(alignment_blocks):
    id2block = dict()

    for blocks in alignment_blocks.values():
        for block in blocks:
            id2block[block.signed_id()] = block

    return id2block


def filter_complement(lengths):
    return [length for i, length in enumerate(sorted(lengths)) if i % 2 == 0]


def get_contig_lengths(assembly_graph):
    return filter_complement([data["length"] for _, _, data in assembly_graph.edges(data=True)])


def get_alignment_block_lengths(alignment_blocks):
    alignment_block_lengths = []

    for blocks in alignment_blocks.values():
        alignment_block_lengths.extend([block.length() for block in blocks])

    return filter_complement(alignment_block_lengths)


def get_path_lengths(paths):
    return [calc_path_length(path[0]) for path in paths]


def calc_path_length(path):
    path_length = 0

    for i in range(len(path)):
        if i % 2 == 0:
            path_length += path[i].length()
        else:
            path_length += sum([seq_length for _, seq_length in path[i]])

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
