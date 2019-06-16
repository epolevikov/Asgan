import os
import argparse

import src.stats as st
import src.paths as ps
import src.hits as ht
import src.aligner as aligner
import src.assembly_graph as asg
import src.adjacency_graph as adg
import src.synteny_blocks as sb
import src.breakpoint_graph as bpg
import src.gfa_parser as gfa_parser
import src.output_generator as out_gen

import networkx as nx


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-query")
    parser.add_argument("--input-target")
    parser.add_argument("--out-dir")
    return parser.parse_args()


def main():
    # Running the pipeline
    args = parse_args()
    os.mkdir(args.out_dir)

    gfa_query, gfa_target = args.input_query, args.input_target

    print("Parsing assembly graphs..")
    assembly_graph_query = asg.parse_assembly_graph(gfa_query)
    assembly_graph_target = asg.parse_assembly_graph(gfa_target)

    repeats_query = asg.get_repeats(assembly_graph_query)
    repeats_target = asg.get_repeats(assembly_graph_target)

    print("Extracting sequences..")
    sequences_query = gfa_parser.extract_sequences(gfa_query, out_dir=args.out_dir,
                                                   out_file="sequences_query.fasta")
    sequences_target = gfa_parser.extract_sequences(gfa_target, out_dir=args.out_dir,
                                                    out_file="sequences_target.fasta")

    print("Aligning sequences..")
    raw_hits = aligner.align(sequences_query, sequences_target, args.out_dir)
    filtered_hits = ht.filter_repeats(raw_hits, repeats_query, repeats_target)
    processed_hits = ht.process_raw_hits(filtered_hits)

    print("Finding shared paths..")
    synteny_blocks_query, synteny_blocks_target = sb.extract_synteny_blocks(processed_hits)

    adjacency_graph_query = adg.build_adjacency_graph(assembly_graph_query, synteny_blocks_query)
    adjacency_graph_target = adg.build_adjacency_graph(assembly_graph_target, synteny_blocks_target)

    breakpoint_graph = bpg.build_breakpoint_graph(adjacency_graph_query, synteny_blocks_query,
                                                  adjacency_graph_target, synteny_blocks_target)

    max_matching = nx.max_weight_matching(breakpoint_graph)

    path_components = bpg.build_path_components(breakpoint_graph, max_matching)
    unused_edges = bpg.get_unused_edges(breakpoint_graph, max_matching)

    number_united_components = bpg.unite_cycles(path_components, unused_edges)

    synteny_paths = ps.build_synteny_paths(path_components)
    path_sequences_query = ps.build_path_sequences(synteny_blocks_query, synteny_paths,
                                                   adjacency_graph_query)
    path_sequences_target = ps.build_path_sequences(synteny_blocks_target, synteny_paths,
                                                    adjacency_graph_target)

    print("Calculating stats..")
    stats = st.calc_stats(assembly_graph_query, synteny_blocks_query, path_sequences_query,
                          assembly_graph_target, synteny_blocks_target, path_sequences_target,
                          synteny_paths, number_united_components, raw_hits, args.out_dir)

    # Generating output
    block_attributes = sb.set_block_attributes(synteny_paths)
    out_gen.adjacency_graph_save_dot(adjacency_graph_query, out_dir=args.out_dir,
                                     block_attributes=block_attributes,
                                     out_file="adjacency_graph_query.gv")
    out_gen.adjacency_graph_save_dot(adjacency_graph_target, out_dir=args.out_dir,
                                     block_attributes=block_attributes,
                                     out_file="adjacency_graph_target.gv")

    out_gen.save_path_sequences(path_sequences_query, path_sequences_target,
                                out_dir=args.out_dir)

    out_gen.output_stats(stats, out_dir=args.out_dir)

    os.remove(sequences_query)
    os.remove(sequences_target)

    '''
    out_gen.save_blocks(synteny_blocks_query, synteny_blocks_target,
                        out_dir=args.out_dir, out_file="synteny_blocks.txt")

    out_gen.assembly_graph_save_dot(graph=assembly_graph_query,
                                    out_file="assembly_graph_query.gv",
                                    out_dir=args.out_dir)
    out_gen.assembly_graph_save_dot(graph=assembly_graph_target,
                                    out_file="assembly_graph_target.gv",
                                    out_dir=args.out_dir)

    out_gen.save_raw_hits(raw_hits, out_dir=args.out_dir, out_file="raw_hits.txt")

    out_gen.breakpoint_graph_save_dot(breakpoint_graph, max_matching,
                                      out_dir=args.out_dir,
                                      out_file="breakpoint_graph.gv")

    out_gen.path_sequences_save_fasta(path_sequences_query, sequences_query,
                                      path_sequences_target, sequences_target,
                                      args.out_dir)
    '''
