import networkx as nx


def calc_stats(assembly_graph_query, assembly_graph_target, aln_blocks_query, aln_blocks_target, paths):
    query_num_wcc = nx.number_weakly_connected_components(assembly_graph_query)
    target_num_wcc = nx.number_weakly_connected_components(assembly_graph_target)

    query_num_contigs = len(assembly_graph_query.edges())
    target_num_contigs = len(assembly_graph_target.edges())

    stats = dict()

    stats["query_num_wcc"] = query_num_wcc
    stats["target_num_wcc"] = target_num_wcc
    stats["query_num_contigs"] = query_num_contigs
    stats["target_num_contigs"] = target_num_contigs
    stats["num_paths"] = len(paths)

    return stats
