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

    def query_hit_length(self):
        return self.query_end - self.query_start

    def target_hit_length(self):
        return self.target_end - self.target_start

    def query_mapping_rate(self):
        return self.query_hit_length() / self.query_len

    def target_mapping_rate(self):
        return self.target_hit_length() / self.target_len

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


def align(contigs_query, contigs_target):
    output_file = "minimap.paf"

    _run_minimap(contigs_query, contigs_target, output_file)

    with open(output_file) as f:
        raw_hits = [RawPafHit(raw_hit) for raw_hit in f]

    os.remove(output_file)
    return raw_hits


def _run_minimap(contigs_query, contigs_target, outfile):
    MINIMAP_BIN = "lib/Flye/bin/flye-minimap2"
    cmd = [MINIMAP_BIN]
    cmd.extend(["--secondary=no"])
    cmd.extend(["-cx", "asm10"])
    cmd.extend([contigs_target, contigs_query])
    cmd.extend([">", outfile])
    cmd.extend(["2> /dev/null"])
    cmd = " ".join(cmd)

    print(cmd)

    os.system(cmd)
