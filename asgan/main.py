import os
import argparse

import asgan.aligner as aligner
import asgan.gfa_parser as gfa_parser
import asgan.output_generator as out_gen
import asgan.hits_processing as hits_proc
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
    hits = hits_proc.process_raw_hits(raw_hits)

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

    aln_blocks_query, aln_blocks_target, num_aln_blocks = align_contigs(args)
    out_gen.output_blocks_info(args.out_dir,
                               aln_blocks_query,
                               aln_blocks_target)

    assembly_graph_query = parse_assembly_graph(
        args.graph_query, aln_blocks_query)
    assembly_graph_target = parse_assembly_graph(
        args.graph_target, aln_blocks_target)

    alignment_graph_query = aln_gr_proc.build_alignment_graph(
        assembly_graph_query, aln_blocks_query)
    alignment_graph_target = aln_gr_proc.build_alignment_graph(
        assembly_graph_target, aln_blocks_target)

    out_gen.alignment_graph_save_dot(alignment_graph_query, args.out_dir,
                                     "alignment_graph_query")
    out_gen.alignment_graph_save_dot(alignment_graph_target, args.out_dir,
                                     "alignment_graph_target")

    breakpoint_graph = bp_gr_proc.build_breakpoint_graph(
        alignment_graph_query, alignment_graph_target, num_aln_blocks)

    max_matching = nx.max_weight_matching(breakpoint_graph)
    out_gen.breakpoint_graph_save_dot(breakpoint_graph, max_matching,
                                      args.out_dir)

    paths = bp_gr_proc.reconstruct_paths(breakpoint_graph, max_matching)

    with open("{}/paths.txt".format(args.out_dir), "w") as f:
        f.write("graph {\n")
        f.write("  edge [penwidth=5]\n")

        for node, data in paths.nodes(data=True):
            f.write("  {} [label=\"{}\"]\n".format(node, data["label"]))

        for (node_from, node_to) in paths.edges():
            f.write("  {} -- {} \n".format(node_from, node_to))

        f.write("}\n")


if __name__ == "__main__":
    main()
