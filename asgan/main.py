import os
import argparse

import asgan.stats as st
import asgan.paths as ps
import asgan.hits as hits
import asgan.aligner as aligner
import asgan.assembly_graph as asg
import asgan.alignment_graph as alg
import asgan.alignment_blocks as alb
import asgan.breakpoint_graph as bpg
import asgan.gfa_parser as gfa_parser
import asgan.output_generator as out_gen

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
    processed_hits = hits.process_raw_hits(raw_hits, args)

    alignment_blocks_query, alignment_blocks_target = alb.extract_alignment_blocks(processed_hits)
    grouped_alignment_blocks_query = alb.group_by_sequence(alignment_blocks_query)
    grouped_alignment_blocks_target = alb.group_by_sequence(alignment_blocks_target)

    return grouped_alignment_blocks_query, grouped_alignment_blocks_target


def parse_assembly_graph(assembly_graph_gfa, alignment_blocks):
    sequences, links = gfa_parser.parse_gfa(assembly_graph_gfa)
    assembly_graph = asg.build_assembly_graph(sequences, links)
    return assembly_graph


def main():
    args = get_args()
    os.mkdir(args.out_dir)

    print("Aligning contigs")
    alignment_blocks_query, alignment_blocks_target = align_contigs(args)
    out_gen.output_blocks_info(alignment_blocks_query, alignment_blocks_target, args.out_dir)

    print("Parsing assembly graphs")
    assembly_graph_query = parse_assembly_graph(args.graph_query, alignment_blocks_query)
    assembly_graph_target = parse_assembly_graph(args.graph_target, alignment_blocks_target)

    out_gen.assembly_graph_save_dot(assembly_graph_query, "assembly_graph_query", args.out_dir)
    out_gen.assembly_graph_save_dot(assembly_graph_target, "assembly_graph_target", args.out_dir)

    print("Building alignment graphs")
    alignment_graph_query = alg.build_alignment_graph(assembly_graph_query, alignment_blocks_query)
    alignment_graph_target = alg.build_alignment_graph(assembly_graph_target, alignment_blocks_target)

    print("Finding shared paths")
    breakpoint_graph = bpg.build_breakpoint_graph(alignment_graph_query, alignment_graph_target)

    max_matching = nx.max_weight_matching(breakpoint_graph)
    out_gen.breakpoint_graph_save_dot(breakpoint_graph, max_matching, args.out_dir)

    path_components, alignment_block_paths = ps.reconstruct_alignment_block_paths(
        alignment_graph_query, alignment_blocks_query,
        alignment_graph_target, alignment_blocks_target,
        breakpoint_graph, max_matching)

    full_paths_query = ps.reconstruct_full_paths(
        alignment_block_paths, alignment_graph_query, alignment_blocks_query)
    full_paths_target = ps.reconstruct_full_paths(
        alignment_block_paths, alignment_graph_target, alignment_blocks_target)

    out_gen.save_full_paths(full_paths_query, full_paths_target, args.out_dir)

    unused_edges = bpg.get_unused_edges(breakpoint_graph, max_matching)
    out_gen.paths_graph_save_dot(path_components, unused_edges, args.out_dir)

    block_colors, block_styles = alb.set_block_attributes(alignment_block_paths)

    out_gen.alignment_graph_save_dot(alignment_graph_query, "alignment_graph_query",
                                     block_colors, block_styles, args.out_dir)
    out_gen.alignment_graph_save_dot(alignment_graph_target, "alignment_graph_target",
                                     block_colors, block_styles, args.out_dir)

    print("Calculating stats")
    stats = st.calc_stats(assembly_graph_query,  alignment_blocks_query, full_paths_query,
                          assembly_graph_target, alignment_blocks_target, full_paths_target,
                          args.out_dir)
    out_gen.output_stats(stats, args.out_dir)


if __name__ == "__main__":
    main()
