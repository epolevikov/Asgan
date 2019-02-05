import sys


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


def main():
    input_file = sys.argv[1]
    output_file = input_file.split(".")[0] + ".gfa"

    with open(sys.argv[1]) as fin, open(output_file, "w") as fout:
        fout.write("H\tVN:Z:1.0\n")
        for header, seq in read_fasta(fin):
            fout.write("S\t{}\t*\tLN:i:{}\n".format(header, len(seq)))


if __name__ == "__main__":
    main()
