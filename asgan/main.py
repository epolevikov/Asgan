import os
import argparse

import asgan.stats as st
import asgan.aligner as aligner
import asgan.gfa_parser as gfa_parser
import asgan.output_generator as out_gen
import asgan.hits_processing as hits_proc
import asgan.paths_processing as paths_proc
import asgan.assembly_graph_processing as agr_proc
import asgan.alignment_blocks_processing as ab_proc
import asgan.alignment_graph_processing as aln_gr_proc
import asgan.breakpoint_graph_processing as bp_gr_proc

import networkx as nx


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph-query")
    parser.add_argument("--graph-target")
    parser.add_argument("--contigs-query")
    parser.add_argument("--contigs-target")
    parser.add_argument("--out-dir")
    return parser.parse_args()


def align_contigs(args):
    raw_hits = aligner.align(args.contigs_query, args.contigs_target)
    hits = hits_proc.process_raw_hits(raw_hits, args)

    aln_blocks_query, aln_blocks_target = \
        ab_proc.extract_alignment_blocks(hits)
    grouped_aln_blocks_query = ab_proc.group_by_sequence(aln_blocks_query)
    grouped_aln_blocks_target = ab_proc.group_by_sequence(aln_blocks_target)

    return (grouped_aln_blocks_query,
            grouped_aln_blocks_target,
            len(aln_blocks_query))


def parse_assembly_graph(graph_gfa, alignment_blocks):
    sequences, links = gfa_parser.parse_gfa(graph_gfa)
    assembly_graph = agr_proc.build_assembly_graph(sequences, links)
    assembly_graph = agr_proc.remove_components_without_alignment_blocks(
                        assembly_graph, alignment_blocks)
    return assembly_graph


def main():
    args = get_args()
    os.mkdir(args.out_dir)

    print("1. Aligning contigs")
    aln_blocks_query, aln_blocks_target, num_aln_blocks = align_contigs(args)
    out_gen.output_blocks_info(args.out_dir,
                               aln_blocks_query,
                               aln_blocks_target)

    print("2. Parsing assembly graphs")
    assembly_graph_query = parse_assembly_graph(
        args.graph_query, aln_blocks_query)
    assembly_graph_target = parse_assembly_graph(
        args.graph_target, aln_blocks_target)

    out_gen.assembly_graph_save_dot(assembly_graph_query, args.out_dir,
                                    "assembly_graph_query")
    out_gen.assembly_graph_save_dot(assembly_graph_target, args.out_dir,
                                    "assembly_graph_target")

    print("3. Building alignment graphs")
    alignment_graph_query = aln_gr_proc.build_alignment_graph(
        assembly_graph_query, aln_blocks_query)
    alignment_graph_target = aln_gr_proc.build_alignment_graph(
        assembly_graph_target, aln_blocks_target)

    print("4. Finding shared paths")
    breakpoint_graph = bp_gr_proc.build_breakpoint_graph(
        alignment_graph_query, alignment_graph_target, num_aln_blocks)

    max_matching = nx.max_weight_matching(breakpoint_graph)
    out_gen.breakpoint_graph_save_dot(breakpoint_graph, max_matching,
                                      args.out_dir)

    paths_graph, paths = paths_proc.reconstruct_paths(breakpoint_graph,
                                                      max_matching)
    unused_edges = bp_gr_proc.get_unused_edges(breakpoint_graph, max_matching)
    out_gen.paths_graph_save_dot(paths_graph, unused_edges, args.out_dir)

    block_colors, block_styles = ab_proc.set_block_attributes(paths)

    out_gen.alignment_graph_save_dot(alignment_graph_query, args.out_dir,
                                     "alignment_graph_query-1", block_colors,
                                     block_styles)
    out_gen.alignment_graph_save_dot(alignment_graph_target, args.out_dir,
                                     "alignment_graph_target-1", block_colors,
                                     block_styles)

    print("5. Calculating stats")
    stats = st.get_stats(assembly_graph_query, aln_blocks_query,
                         assembly_graph_target, aln_blocks_target,
                         paths)
    out_gen.output_stats(stats, args.out_dir)

    # out_gen.alignment_graphs_save_dot(alignment_graph_query,
    #                                  alignment_graph_target, block_colors,
    #                                  args.out_dir, "alignment_graphs.gv")


if __name__ == "__main__":
    main()
