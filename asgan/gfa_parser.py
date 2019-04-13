from collections import namedtuple

Link = namedtuple("Link", ["from_name", "from_strand", "to_name", "to_strand"])
Sequence = namedtuple("Sequence", ["name", "length", "seq", "is_repeat"])


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


def extract_sequences(gfa_file, out_dir, out_file):
    out_file = "{}/{}".format(out_dir, out_file)

    with open(gfa_file) as fin, open(out_file, "w") as fout:
        for line in fin:
            record = line.strip().split()
            record_type = record[0]

            if record_type == RecordType.SEQUENCE:
                sequence = _parse_sequence(record)
                fout.write(">{}\n{}\n".format(sequence.name, sequence.seq))

    return out_file


def _parse_sequence(record):
    name = record[1]
    seq = record[2]
    length = len(seq)
    is_repeat = False

    if len(record) > 3 and record[3].startswith("r"):
        is_repeat = (int(record[3].split(":")[-1]) == 1)

    return Sequence(name, length, seq, is_repeat)


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
