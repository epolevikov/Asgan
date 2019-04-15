import networkx as nx


def calc_stats(assembly_graph_query, synteny_blocks_query, path_sequences_query,
               assembly_graph_target, synteny_blocks_target, path_sequences_target,
               synteny_paths, number_united_components, out_dir):
    number_wcc_query = nx.number_weakly_connected_components(assembly_graph_query)
    number_wcc_target = nx.number_weakly_connected_components(assembly_graph_target)

    # sequences
    sequence_lengths_query = calc_sequence_lengths(assembly_graph_query)
    sequence_lengths_target = calc_sequence_lengths(assembly_graph_target)

    number_sequences_query = len(sequence_lengths_query)
    number_sequences_target = len(sequence_lengths_target)

    sequences_total_length_query = sum(sequence_lengths_query)
    sequences_total_length_target = sum(sequence_lengths_target)

    sequences_n50_query, sequences_l50_query = calc_nx(sequence_lengths_query)
    sequences_n50_target, sequences_l50_target = calc_nx(sequence_lengths_target)

    # synteny blocks
    block_lengths_query = calc_synteny_block_lengths(synteny_blocks_query)
    block_lengths_target = calc_synteny_block_lengths(synteny_blocks_target)

    blocks_total_length_query = sum(block_lengths_query)
    blocks_total_length_target = sum(block_lengths_target)

    number_blocks = len(block_lengths_query)

    blocks_n50_query, blocks_l50_query = calc_nx(block_lengths_query)
    blocks_n50_target, blocks_l50_target = calc_nx(block_lengths_target)

    # synteny paths
    path_lengths_query = calc_path_lengths(path_sequences_query)
    path_lengths_target = calc_path_lengths(path_sequences_target)

    paths_total_length_query = sum(path_lengths_query)
    paths_total_length_target = sum(path_lengths_target)

    number_paths = len(path_lengths_query)

    paths_n50_query, paths_l50_query = calc_nx(path_lengths_query)
    paths_n50_target, paths_l50_target = calc_nx(path_lengths_target)

    # link types:
    link_types = calc_link_types(synteny_paths, synteny_blocks_query, synteny_blocks_target)

    stats = dict()

    stats["link_types"] = link_types

    stats["number_wcc_query"] = number_wcc_query
    stats["number_wcc_target"] = number_wcc_target

    stats["number_sequences_query"] = number_sequences_query
    stats["number_sequences_target"] = number_sequences_target
    stats["sequences_n50_query"] = sequences_n50_query
    stats["sequences_n50_target"] = sequences_n50_target
    stats["sequences_l50_query"] = sequences_l50_query
    stats["sequences_l50_target"] = sequences_l50_target
    stats["sequences_total_length_query"] = sequences_total_length_query
    stats["sequences_total_length_target"] = sequences_total_length_target

    stats["number_blocks"] = number_blocks
    stats["blocks_n50_query"] = blocks_n50_query
    stats["blocks_n50_target"] = blocks_n50_target
    stats["blocks_l50_query"] = blocks_l50_query
    stats["blocks_l50_target"] = blocks_l50_target
    stats["blocks_total_length_query"] = blocks_total_length_query
    stats["blocks_total_length_target"] = blocks_total_length_target

    stats["number_paths"] = number_paths
    stats["paths_n50_query"] = paths_n50_query
    stats["paths_n50_target"] = paths_n50_target
    stats["paths_l50_query"] = paths_l50_query
    stats["paths_l50_target"] = paths_l50_target
    stats["paths_total_length_query"] = paths_total_length_query
    stats["paths_total_length_target"] = paths_total_length_target

    stats["link_types"] = link_types
    stats["number_united_components"] = number_united_components

    return stats


def calc_link_types(synteny_paths, synteny_blocks_query, synteny_blocks_target):
    def get_link_type(block_from, block_to, id2block_query, id2block_target):
        link_types = [[0, 2], [1, 3]]

        seq_from_query = id2block_query[block_from].sequence_name
        seq_to_query = id2block_query[block_to].sequence_name

        seq_from_target = id2block_target[block_from].sequence_name
        seq_to_target = id2block_target[block_to].sequence_name

        row = seq_from_query != seq_to_query
        col = seq_from_target != seq_to_target

        return link_types[row][col]

    id2block_query = build_id2block_dict(synteny_blocks_query)
    id2block_target = build_id2block_dict(synteny_blocks_target)
    link_types = [0, 0, 0, 0]

    for path in synteny_paths:
        if len(path) > 1:
            for i in range(len(path) - 1):
                block_from, block_to = path[i], path[i + 1]
                link_type = get_link_type(block_from, block_to,
                                          id2block_query, id2block_target)
                link_types[link_type] += 1

    return link_types


def build_id2block_dict(synteny_blocks):
    id2block = dict()

    for blocks in synteny_blocks.values():
        for block in blocks:
            id2block[block.signed_id()] = block

    return id2block


def filter_complement(lengths):
    return [length for i, length in enumerate(sorted(lengths)) if i % 2 == 0]


def calc_sequence_lengths(assembly_graph):
    return filter_complement([data["length"] for _, _, data in assembly_graph.edges(data=True)])


def calc_synteny_block_lengths(synteny_blocks):
    synteny_block_lengths = []

    for blocks in synteny_blocks.values():
        synteny_block_lengths.extend([block.length() for block in blocks])

    return filter_complement(synteny_block_lengths)


def calc_path_lengths(paths):
    return [calc_path_length(path) for path in paths]


def calc_path_length(path):
    path_length = 0

    for subpath in path:
        if isinstance(subpath, list):
            path_length += sum([block.length() for block in subpath])
        else:
            path_length += subpath.length()

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
