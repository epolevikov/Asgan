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
    parser.add_argument("--out-file")
    args = parser.parse_args()

    with open(args.input_fasta) as fin, open(args.out_file, "w") as fout:
        fout.write("H\tVN:Z:1.0\n")
        for header, seq in read_fasta(fin):
            fout.write("S\t{}\t{}\tLN:i:{}\n".format(header, seq, len(seq)))


if __name__ == "__main__":
    main()
