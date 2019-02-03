import networkx as nx


class Stats:
    def __init__(self, query_num_wcc, query_num_contigs,
                 target_num_wcc, target_num_contigs, paths):
        self.query_num_wcc = query_num_wcc
        self.query_num_contigs = query_num_contigs
        self.target_num_wcc = target_num_wcc
        self.target_num_contigs = target_num_contigs
        self.paths = paths


def get_stats(graph_query, graph_target, paths):
    query_num_wcc = nx.number_weakly_connected_components(graph_query)
    query_num_contigs = len(graph_query.edges())
    target_num_wcc = nx.number_weakly_connected_components(graph_target)
    target_num_contigs = len(graph_target.edges())

    return Stats(query_num_wcc, query_num_contigs,
                 target_num_wcc, target_num_contigs,
                 paths)
