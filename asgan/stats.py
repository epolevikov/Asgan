import networkx as nx


def calc_stats(assembly_graph_query, synteny_blocks_query, path_sequences_query,
               assembly_graph_target, synteny_blocks_target, path_sequences_target,
               synteny_paths, number_united_components, raw_hits, out_dir):
    stats = dict()

    # number wcc
    number_wcc_query = number_wcc(assembly_graph_query, synteny_blocks_query)
    number_wcc_target = number_wcc(assembly_graph_target, synteny_blocks_target)

    stats["number_wcc_query"] = number_wcc_query
    stats["number_wcc_target"] = number_wcc_target

    # sequences
    sequence_lengths_query = calc_sequence_lengths(assembly_graph_query,
                                                   synteny_blocks_query)
    sequence_lengths_target = calc_sequence_lengths(assembly_graph_target,
                                                    synteny_blocks_target)

    number_sequences_query = len(sequence_lengths_query)
    number_sequences_target = len(sequence_lengths_target)

    sequences_total_length_query = sum(sequence_lengths_query)
    sequences_total_length_target = sum(sequence_lengths_target)

    sequences_n50_query, sequences_l50_query = calc_nx(sequence_lengths_query)
    sequences_n50_target, sequences_l50_target = calc_nx(sequence_lengths_target)

    stats["number_sequences_query"] = number_sequences_query
    stats["number_sequences_target"] = number_sequences_target
    stats["sequences_n50_query"] = sequences_n50_query
    stats["sequences_n50_target"] = sequences_n50_target
    stats["sequences_l50_query"] = sequences_l50_query
    stats["sequences_l50_target"] = sequences_l50_target
    stats["sequences_total_length_query"] = sequences_total_length_query
    stats["sequences_total_length_target"] = sequences_total_length_target

    # synteny blocks
    block_lengths_query = calc_synteny_block_lengths(synteny_blocks_query)
    block_lengths_target = calc_synteny_block_lengths(synteny_blocks_target)

    blocks_total_length_query = sum(block_lengths_query)
    blocks_total_length_target = sum(block_lengths_target)

    number_blocks = len(block_lengths_query)

    blocks_n50_query, blocks_l50_query = calc_nx(block_lengths_query)
    blocks_n50_target, blocks_l50_target = calc_nx(block_lengths_target)

    stats["number_blocks"] = number_blocks
    stats["blocks_n50_query"] = blocks_n50_query
    stats["blocks_n50_target"] = blocks_n50_target
    stats["blocks_l50_query"] = blocks_l50_query
    stats["blocks_l50_target"] = blocks_l50_target
    stats["blocks_total_length_query"] = blocks_total_length_query
    stats["blocks_total_length_target"] = blocks_total_length_target

    # synteny paths
    path_lengths_query = calc_path_lengths(path_sequences_query)
    path_lengths_target = calc_path_lengths(path_sequences_target)

    paths_total_length_query = sum(path_lengths_query)
    paths_total_length_target = sum(path_lengths_target)

    number_paths = len(path_lengths_query)

    paths_n50_query, paths_l50_query = calc_nx(path_lengths_query)
    paths_n50_target, paths_l50_target = calc_nx(path_lengths_target)

    stats["number_paths"] = number_paths
    stats["paths_n50_query"] = paths_n50_query
    stats["paths_n50_target"] = paths_n50_target
    stats["paths_l50_query"] = paths_l50_query
    stats["paths_l50_target"] = paths_l50_target
    stats["paths_total_length_query"] = paths_total_length_query
    stats["paths_total_length_target"] = paths_total_length_target

    # link types:
    link_types = calc_link_types(synteny_paths, synteny_blocks_query, synteny_blocks_target)
    stats["link_types"] = link_types

    # united components:
    stats["number_united_components"] = number_united_components

    # alignment identity

    mean_alignment_identity = calc_mean_alignment_identity(raw_hits)
    total_alignment_identity = calc_total_alignment_identity(raw_hits)

    stats["mean_alignment_identity"] = mean_alignment_identity
    stats["total_alignment_identity"] = total_alignment_identity

    # assembly coverage

    query_assembly_coverage, target_assembly_coverage = calc_assembly_coverage(
        raw_hits, assembly_graph_query, assembly_graph_target)

    stats["query_assembly_coverage"] = query_assembly_coverage
    stats["target_assembly_coverage"] = target_assembly_coverage

    return stats


def number_wcc(assembly_graph, synteny_blocks):
    number_wcc = 0

    for component in nx.weakly_connected_component_subgraphs(assembly_graph, copy=False):
        if contains_synteny_blocks(component, synteny_blocks):
            number_wcc += 1

    return number_wcc


def calc_sequence_lengths(assembly_graph, synteny_blocks):
    sequence_lengths = []

    for component in nx.weakly_connected_component_subgraphs(assembly_graph, copy=False):
        if contains_synteny_blocks(component, synteny_blocks):
            sequence_lengths.extend([data["length"] for (_, _, data)
                                     in component.edges(data=True)])

    return filter_complement(sequence_lengths)


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


def calc_mean_alignment_identity(raw_hits):
    alignment_identities = [hit.alignment_identity() for hit in raw_hits]
    mean_alignment_identity = sum(alignment_identities) / len(alignment_identities)
    return round(mean_alignment_identity, 3)


def calc_total_alignment_identity(raw_hits):
    matching_bases = [hit.matching_bases for hit in raw_hits]
    number_bases = [hit.number_bases for hit in raw_hits]
    total_alignment_identity = sum(matching_bases) / sum(number_bases)
    return round(total_alignment_identity, 3)


def calc_assembly_coverage(raw_hits, assembly_graph_query, assembly_graph_target):
    query_hit_lengths = [hit.query_hit_length() for hit in raw_hits]
    target_hit_lengths = [hit.target_hit_length() for hit in raw_hits]

    query_sequence_lengths = [data["length"] for (_, _, data)
                              in assembly_graph_query.edges(data=True)]
    target_sequence_lengths = [data["length"] for (_, _, data)
                               in assembly_graph_target.edges(data=True)]

    query_hit_total_length = sum(query_hit_lengths)
    query_sequence_total_length = sum(query_sequence_lengths) / 2
    query_assembly_coverage = query_hit_total_length / query_sequence_total_length
    query_assembly_coverage = round(query_assembly_coverage, 3)

    target_hit_total_length = sum(target_hit_lengths)
    target_sequence_total_length = sum(target_sequence_lengths) / 2
    target_assembly_coverage = target_hit_total_length / target_sequence_total_length
    target_assembly_coverage = round(target_assembly_coverage, 3)

    return query_assembly_coverage, target_assembly_coverage


def filter_complement(lengths):
    return [length for i, length in enumerate(sorted(lengths)) if i % 2 == 0]


def contains_synteny_blocks(component, synteny_blocks):
    for (_, _, data) in component.edges(data=True):
        if data["name"] in synteny_blocks:
            return True

    return False


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
