import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-gfa")
    parser.add_argument("--out-file")
    args = parser.parse_args()

    with open(args.input_gfa) as fin, open(args.out_file, "w") as fout:
        for line in fin:
            record = line.strip().split()
            if record[0] == "S":
                fout.write(">{}\n{}\n".format(record[1], record[2]))


if __name__ == "__main__":
    main()
