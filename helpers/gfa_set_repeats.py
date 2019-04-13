import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-gfa")
    parser.add_argument("--out-file")
    args = parser.parse_args()

    depths = []
    lengths = []

    with open(args.input_gfa) as f:
        for line in f:
            record = line.strip().split()
            if record[0] == "S":
                depth = int(record[3].split(":")[-1])

                depths.append(depth)
                lengths.append(len(record[2]))

    depths_weighted_sum = sum([depths[i] * lengths[i] for i in range(len(depths))])
    depths_weighted_mean = depths_weighted_sum / sum(lengths)

    with open(args.input_gfa) as fin, open(args.out_file, "w") as fout:
        for line in fin:
            record = line.strip().split()
            if record[0] == "S":
                depth = int(record[3].split(":")[-1])
                is_repeat = (depth > 1.75 * depths_weighted_mean)
                fout.write("S\t{}\t{}\tr:i:{}\n".format(record[1], record[2], int(is_repeat)))
            else:
                fout.write(line)


if __name__ == "__main__":
    main()
