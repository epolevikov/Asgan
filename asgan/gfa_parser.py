

class RecordType:
    SEQUENCE = "S"
    LINK = "L"


class Link:
    def __init__(self, from_name, from_strand, to_name, to_strand):
        self.from_name = from_name
        self.from_strand = from_strand
        self.to_name = to_name
        self.to_strand = to_strand


class Sequence:
    def __init__(self, name, length):
        self.name = name
        self.length = length


def parse_gfa(gfa_file):
    sequences, links = [], []

    with open(gfa_file) as f:
        for line in f:
            record = line.strip().split()
            record_type = record[0]

            if record_type == RecordType.SEQUENCE:
                sequences.append(_parse_sequence(record))

            if record_type == RecordType.LINK:
                links.append(_parse_link(record))

    return sequences, links


def _parse_sequence(record):
    name, seq = record[1], record[2]
    length = record[3].split(":")[-1] if seq == "*" else len(seq)
    return Sequence(name, int(length))


def _parse_link(record):
    from_name, from_strand, to_name, to_strand = record[1:5]
    return Link(from_name, from_strand, to_name, to_strand)
