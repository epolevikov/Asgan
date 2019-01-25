import os


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

    def __str__(self):
        query_info = "{}\t{}\t{}\t{}\t".format(
            self.query_name, self.query_len,
            self.query_start, self.query_end)

        target_info = "\t{}\t{}\t{}\t{}".format(
            self.target_name, self.target_len,
            self.target_start, self.target_end)

        strand = self.strand
        return query_info + strand + target_info


def align(contigs_query, contigs_target):
    output_file = "minimap.paf"
    _run_minimap(contigs_query, contigs_target, output_file)

    with open(output_file) as f:
        raw_hits = [RawPafHit(raw_hit) for raw_hit in f]

    os.remove(output_file)
    return raw_hits


def _run_minimap(contigs_query, contigs_target, output_file):
    MINIMAP_BIN = "lib/minimap2/minimap2"
    cmd = [MINIMAP_BIN]
    cmd.extend(["-secondary=no"])
    cmd.extend(["-cx", "asm10"])
    cmd.extend([contigs_target, contigs_query])
    cmd.extend([">", output_file])
    cmd.extend(["2> /dev/null"])
    os.system(" ".join(cmd))
