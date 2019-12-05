import os
from asgan.output_generator import pretty_number


class RawPafHit:
    def __init__(self, raw_hit):
        raw_hit = raw_hit.strip().split()

        self.query_name = raw_hit[0]
        self.query_len = int(raw_hit[1])
        self.query_start = int(raw_hit[2])
        self.query_end = int(raw_hit[3])

        self.strand = raw_hit[4]

        self.target_name = raw_hit[5]
        self.target_len = int(raw_hit[6])
        self.target_start = int(raw_hit[7])
        self.target_end = int(raw_hit[8])

        self.matching_bases = int(raw_hit[9])
        self.number_bases = int(raw_hit[10])

    def query_hit_length(self):
        return self.query_end - self.query_start

    def target_hit_length(self):
        return self.target_end - self.target_start

    def alignment_identity(self):
        return float(self.matching_bases) / float(self.number_bases)

    def __str__(self):
        query_info = "{}\t{}\t{}\t{}".format(
            self.query_name,
            pretty_number(self.query_len),
            pretty_number(self.query_start),
            pretty_number(self.query_end))

        target_info = "{}\t{}\t{}\t{}".format(
            self.target_name,
            pretty_number(self.target_len),
            pretty_number(self.target_start),
            pretty_number(self.target_end))

        return "{}\t{}\t{}".format(query_info, self.strand, target_info)


def align(sequences_query, sequences_target, args):
    out_file = args.out_dir + "/minimap.paf"

    run_minimap(sequences_query, sequences_target, out_file, args)

    with open(out_file) as f:
        raw_hits = [RawPafHit(raw_hit) for raw_hit in f]

    os.remove(out_file)
    return raw_hits


def run_minimap(sequences_query, sequences_target, out_file, args):
    asgan_root = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    MINIMAP_BIN = os.path.join(asgan_root, "lib/minimap2/minimap2")

    cmd = [MINIMAP_BIN]
    cmd.extend(["--secondary=no"])
    cmd.extend(["-cx", args.minimap_preset])
    cmd.extend([sequences_target, sequences_query])
    cmd.extend([">", out_file])
    cmd.extend(["2> /dev/null"])
    cmd = " ".join(cmd)

    os.system(cmd)
