import os
import argparse

import asgan.stats as st
import asgan.paths as ps
import asgan.hits as hits
import asgan.flye_repeat as fr
import asgan.aligner as aligner
import asgan.assembly_graph as asg
import asgan.alignment_graph as alg
import asgan.alignment_blocks as alb
import asgan.breakpoint_graph as bpg
import asgan.gfa_parser as gfa_parser
import asgan.output_generator as out_gen

import networkx as nx


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-query")
    parser.add_argument("--input-target")
    parser.add_argument("--ref", type=bool, default=False)
    parser.add_argument("--single-graph", type=bool, default=False)
    parser.add_argument("--out-dir")
    return parser.parse_args()


def parse_input(args):
    gfa_query = args.input_query
    gfa_target = args.input_target

    if args.input_query.endswith(".fasta"):
        gfa_query = fr.run_flye_repeat(
            sequences_fasta=args.input_query,
            out_dir=args.out_dir,
            out_name="graph_query.gfa")

    if args.input_target.endswith(".fasta"):
        gfa_target = fr.run_flye_repeat(
            sequences_fasta=args.input_target,
            out_dir=args.out_dir,
            out_name="graph_target.gfa")

    return gfa_query, gfa_target


def build_assembly_graph_from_gfa(gfa_file):
    sequences, links = gfa_parser.parse_gfa(gfa_file)
    return asg.build_assembly_graph(sequences, links)


def extract_alignment_blocks_from_hits(hits):
    alignment_blocks_query, alignment_blocks_target = alb.extract_alignment_blocks(hits)
    grouped_alignment_blocks_query = alb.group_by_sequence(alignment_blocks_query)
    grouped_alignment_blocks_target = alb.group_by_sequence(alignment_blocks_target)

    return grouped_alignment_blocks_query, grouped_alignment_blocks_target


def main():
    args = parse_args()
    os.mkdir(args.out_dir)

    ''' The code below needs to be refactored since it is not the best way
        to process different input modes. But for now, it's ok. It will be
        improved further.
    '''

    if args.single_graph:
        gfa_query = args.input_query
        assembly_graph_query = build_assembly_graph_from_gfa(gfa_query)
        asg.mark_repeats(assembly_graph_query, normalize_depth=True)
        repeats_query = asg.get_repeats(assembly_graph_query)

        out_gen.assembly_graph_save_dot(assembly_graph_query, "assembly_graph_query.gv", args.out_dir)

        alignment_blocks_query = alb.build_from_sequences(assembly_graph_query, repeats_query)
        alignment_graph_query = alg.build_alignment_graph(assembly_graph_query, alignment_blocks_query)

        breakpoint_graph = bpg.build_breakpoint_graph(alignment_graph_query, alignment_blocks_query,
                                                      alignment_graph_query, alignment_blocks_query)

        max_matching = nx.max_weight_matching(breakpoint_graph)
        out_gen.breakpoint_graph_save_dot(breakpoint_graph, max_matching, args.out_dir)

        path_components = bpg.build_path_components(breakpoint_graph, max_matching)
        unused_edges = bpg.get_unused_edges(breakpoint_graph, max_matching)

        bpg.unite_cycles(path_components, unused_edges, args)

        alignment_block_paths = ps.reconstruct_alignment_block_paths(
            alignment_graph_query, alignment_blocks_query,
            alignment_graph_query, alignment_blocks_query,
            path_components)

        full_paths_query = ps.reconstruct_full_paths(
            alignment_block_paths, alignment_graph_query, alignment_blocks_query)

        out_gen.save_full_paths(full_paths_query, full_paths_query, args.out_dir)
        out_gen.paths_graph_save_dot(path_components, unused_edges, args.out_dir)

        block_colors, block_styles = alb.set_block_attributes(alignment_block_paths)

        out_gen.alignment_graph_save_dot(alignment_graph_query, "alignment_graph_query.gv",
                                         block_colors, block_styles, args.out_dir)

        # print("Calculating stats")
        stats = st.calc_stats(assembly_graph_query,  alignment_blocks_query, full_paths_query,
                              assembly_graph_query, alignment_blocks_query, full_paths_query,
                              args.out_dir)

        out_gen.output_stats(stats, args.out_dir)

        return
    elif args.ref:
        #gfa_query = args.input_query
        gfa_query = gfa_parser.build_gfa_from_fasta(
            sequences_fasta=args.input_query,
            out_dir=args.out_dir,
            out_name="graph_query.gfa")
        gfa_target = gfa_parser.build_gfa_from_fasta(
            sequences_fasta=args.input_target,
            out_dir=args.out_dir,
            out_name="graph_target.gfa")

        print("Building assembly graphs")

        assembly_graph_query = build_assembly_graph_from_gfa(gfa_query)
        assembly_graph_target = build_assembly_graph_from_gfa(gfa_target)

        #asg.mark_repeats(assembly_graph_query, normalize_depth=True)

        #repeats_query = asg.get_repeats(assembly_graph_query)
        repeats_query = set()
        repeats_target = set()

        #sequences_fasta_query = gfa_parser.extract_sequences_from_gfa(
        #    gfa_file=gfa_query, out_dir=args.out_dir,
        #    out_name="sequences_query.fasta")

        sequences_fasta_query = args.input_query
        sequences_fasta_target = args.input_target

        #if args.input_query.endswith(".fasta"):
        #    os.remove(gfa_query)

        # os.remove(gfa_target)
    else:
        gfa_query, gfa_target = parse_input(args)

        print("Building assembly graphs")

        assembly_graph_query = build_assembly_graph_from_gfa(gfa_query)
        assembly_graph_target = build_assembly_graph_from_gfa(gfa_target)

        asg.mark_repeats(assembly_graph_query, normalize_depth=True)
        asg.mark_repeats(assembly_graph_target, normalize_depth=True)

        repeats_query = asg.get_repeats(assembly_graph_query)
        repeats_target = asg.get_repeats(assembly_graph_target)

        sequences_fasta_query = gfa_parser.extract_sequences_from_gfa(
            gfa_file=gfa_query, out_dir=args.out_dir,
            out_name="sequences_query.fasta")

        sequences_fasta_target = gfa_parser.extract_sequences_from_gfa(
            gfa_file=gfa_target, out_dir=args.out_dir,
            out_name="sequences_target.fasta")

        #if args.input_query.endswith(".fasta"):
        #    os.remove(gfa_query)

        #if args.input_target.endswith(".fasta"):
        #    os.remove(gfa_target)

    print("Aligning sequences")

    raw_hits = aligner.align(sequences_fasta_query, sequences_fasta_target)

    raw_hits.sort(key=lambda hit: (hit.query_name, hit.query_start))

    with open("{}/raw_hits.txt".format(args.out_dir), "w") as f:
        for raw_hit in raw_hits:
            f.write(str(raw_hit) + "\n")

    processed_hits = hits.process_raw_hits(raw_hits, repeats_query, repeats_target, args)

    alignment_blocks_query, alignment_blocks_target = extract_alignment_blocks_from_hits(processed_hits)

    asg.mark_sequences_without_alignment_blocks(assembly_graph_query, alignment_blocks_query)
    asg.mark_sequences_without_alignment_blocks(assembly_graph_target, alignment_blocks_target)

    out_gen.assembly_graph_save_dot(assembly_graph_query, "assembly_graph_query.gv", args.out_dir)
    out_gen.assembly_graph_save_dot(assembly_graph_target, "assembly_graph_target.gv", args.out_dir)

    out_gen.output_blocks_info(alignment_blocks_query, alignment_blocks_target, args.out_dir)

    # os.remove(sequences_fasta_query)
    # os.remove(sequences_fasta_target)

    print("Building alignment graphs")

    alignment_graph_query = alg.build_alignment_graph(assembly_graph_query, alignment_blocks_query)
    alignment_graph_target = alg.build_alignment_graph(assembly_graph_target, alignment_blocks_target)

    print("Finding shared paths")

    breakpoint_graph = bpg.build_breakpoint_graph(alignment_graph_query, alignment_blocks_query,
                                                  alignment_graph_target, alignment_blocks_target)

    max_matching = nx.max_weight_matching(breakpoint_graph)

    out_gen.breakpoint_graph_save_dot(breakpoint_graph, max_matching, args.out_dir)

    path_components = bpg.build_path_components(breakpoint_graph, max_matching)
    unused_edges = bpg.get_unused_edges(breakpoint_graph, max_matching)

    bpg.unite_cycles(path_components, unused_edges, args)

    alignment_block_paths = ps.reconstruct_alignment_block_paths(
        alignment_graph_query, alignment_blocks_query,
        alignment_graph_target, alignment_blocks_target,
        path_components)

    full_paths_query = ps.reconstruct_full_paths(
        alignment_block_paths, alignment_graph_query, alignment_blocks_query)
    full_paths_target = ps.reconstruct_full_paths(
        alignment_block_paths, alignment_graph_target, alignment_blocks_target)

    out_gen.save_full_paths(full_paths_query, full_paths_target, args.out_dir)
    out_gen.paths_graph_save_dot(path_components, unused_edges, args.out_dir)

    block_colors, block_styles = alb.set_block_attributes(alignment_block_paths)

    out_gen.alignment_graph_save_dot(alignment_graph_query, "alignment_graph_query.gv",
                                     block_colors, block_styles, args.out_dir)
    out_gen.alignment_graph_save_dot(alignment_graph_target, "alignment_graph_target.gv",
                                     block_colors, block_styles, args.out_dir)

    print("Calculating stats")
    stats = st.calc_stats(assembly_graph_query,  alignment_blocks_query, full_paths_query,
                          assembly_graph_target, alignment_blocks_target, full_paths_target,
                          args.out_dir)

    out_gen.output_stats(stats, args.out_dir)


if __name__ == "__main__":
    main()
