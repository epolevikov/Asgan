

def read_fasta(file_handle):
    header = None
    seq = []

    for line in file_handle:
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


def make_fasta_dict(input_fasta):
    sequences = dict()

    with open(input_fasta) as f:
        for header, seq in read_fasta(f):
            sequences[header] = seq

    return sequences
