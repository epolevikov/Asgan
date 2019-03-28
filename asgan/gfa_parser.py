from collections import namedtuple

Link = namedtuple("Link", ["from_name", "from_strand", "to_name", "to_strand"])
Sequence = namedtuple("Sequence", ["name", "length", "seq", "depth"])


class RecordType:
    SEQUENCE = "S"
    LINK = "L"


def parse_gfa(gfa_file):
    sequences, links = [], []

    with open(gfa_file) as f:
        for line in f:
            record = line.strip().split()
            record_type = record[0]

            if record_type == RecordType.SEQUENCE:
                sequences.append(_parse_sequence(record))

            if record_type == RecordType.LINK:
                link = _parse_link(record)
                links.append(link)
                links.append(_inv_link(link))

    return sequences, list(set(links))


def extract_sequences_from_gfa(gfa_file, out_dir, out_name):
    outfile = "{}/{}".format(out_dir, out_name)

    with open(gfa_file) as fin, open(outfile, "w") as fout:
        for line in fin:
            record = line.strip().split()
            record_type = record[0]

            if record_type == RecordType.SEQUENCE:
                sequence = _parse_sequence(record)
                fout.write(">{}\n{}\n".format(sequence.name, sequence.seq))

    return outfile


def build_gfa_from_fasta(sequences_fasta, out_dir, out_name):
    outfile = "{}/{}".format(out_dir, out_name)

    with open(sequences_fasta) as fin, open(outfile, "w") as fout:
        fout.write("H\tVN:Z:1.0\n")
        for header, seq in read_fasta(fin):
            fout.write("S\t{}\t*\tLN:i:{}\n".format(header, len(seq)))

    return outfile


def read_fasta(filename):
    header = None
    seq = []

    for line in filename:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if header:
                yield header, "".join(seq)
                seq = []
            header = line[1:].split(" ")[0]
        else:
            seq.append(line)

    if header and len(seq):
        yield header, "".join(seq)


def _parse_sequence(record):
    name = record[1]
    seq = record[2]

    if seq == "*":
        length = int(record[3].split(":")[-1])
        depth = 0
    else:
        length = len(seq)
        depth = int(record[3].split(":")[-1])

    return Sequence(name, length, seq, depth)


def _parse_link(record):
    from_name, from_strand, to_name, to_strand = record[1:5]
    return Link(from_name, from_strand, to_name, to_strand)


def _inv_link(link):
    def inv_strand(strand):
        return ["+", "-"][strand == "+"]

    from_name = link.to_name
    from_strand = inv_strand(link.to_strand)
    to_name = link.from_name
    to_strand = inv_strand(link.from_strand)

    return Link(from_name, from_strand, to_name, to_strand)
