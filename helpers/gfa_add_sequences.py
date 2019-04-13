import argparse


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta")
    parser.add_argument("--input-gfa")
    parser.add_argument("--out-file")
    args = parser.parse_args()

    sequences = dict()

    with open(args.input_fasta) as f:
        for header, seq in read_fasta(f):
            sequences[header] = seq

    with open(args.input_gfa) as fin, open(args.out_file, "w") as fout:
        for line in fin:
            record = line.strip().split()
            if record[0] == "S":
                sequence = sequences[record[1]]
                fout.write("S\t{}\t{}\tLN:i{}\n".format(record[1], sequence, len(sequence)))
            else:
                fout.write(line)


if __name__ == "__main__":
    main()
