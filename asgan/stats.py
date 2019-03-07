import networkx as nx


class Stats:
    def __init__(self, query_num_wcc, query_num_contigs, query_paths_coverage,
                 target_num_wcc, target_num_contigs, target_paths_coverage,
                 paths):
        self.query_num_wcc = query_num_wcc
        self.query_num_contigs = query_num_contigs
        self.query_paths_coverage = query_paths_coverage

        self.target_num_wcc = target_num_wcc
        self.target_num_contigs = target_num_contigs
        self.target_paths_coverage = target_paths_coverage

        self.paths = paths


def get_stats(graph_query, aln_blocks_query, graph_target, aln_blocks_target,
              paths):
    query_num_wcc = nx.number_weakly_connected_components(graph_query)
    query_num_contigs = len(graph_query.edges())
    query_paths_coverage = calculate_paths_coverage(aln_blocks_query)

    target_num_wcc = nx.number_weakly_connected_components(graph_target)
    target_num_contigs = len(graph_target.edges())
    target_paths_coverage = calculate_paths_coverage(aln_blocks_target)

    return Stats(query_num_wcc, query_num_contigs, query_paths_coverage,
                 target_num_wcc, target_num_contigs, target_paths_coverage,
                 paths)


def calculate_paths_coverage(alignment_blocks):
    contigs_total_len = 0
    blocks_total_len = 0

    for contig, blocks in alignment_blocks.items():
        contigs_total_len += blocks[0].seq_length
        for block in blocks:
            blocks_total_len += block.length()

    return round(blocks_total_len / contigs_total_len, 2)
