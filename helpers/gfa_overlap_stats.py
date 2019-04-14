import argparse


def cigar_total_length(cigar_string, match_only=False):
    total_length = 0
    curr_length = ""

    for ch in cigar_string:
        if ch.isdigit():
            curr_length += ch
        else:
            if not match_only or (match_only and ch == "M"):
                total_length += int(curr_length)

            curr_length = ""

    return total_length


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-gfa")
    args = parser.parse_args()

    overlap_lengths = []
    identities = []

    with open(args.input_gfa) as f:
        for line in f:
            line = line.strip().split()

            if line[0] == "L":
                total_length = cigar_total_length(line[5])
                match_length = cigar_total_length(line[5], match_only=True)
                overlap_lengths.append(total_length)
                identities.append(match_length / total_length)

    overlap_lengths.sort()
    min_overlap = overlap_lengths[0]
    max_overlap = overlap_lengths[-1]
    median_overlap = overlap_lengths[len(overlap_lengths) // 2]
    mean_overlap = sum(overlap_lengths) // len(overlap_lengths)

    print("Overlap stats:")
    print("min: {}".format(min_overlap))
    print("max: {}".format(max_overlap))
    print("median: {}".format(median_overlap))
    print("mean: {}\n".format(mean_overlap))

    identities.sort()
    min_identity = round(identities[0], 3)
    max_identity = round(identities[-1], 3)
    median_identity = round(identities[len(identities) // 2], 3)
    mean_identity = round(sum(identities) / len(identities), 3)

    print("Identity stats:")
    print("min: {}".format(min_identity))
    print("max: {}".format(max_identity))
    print("median: {}".format(median_identity))
    print("mean: {}".format(mean_identity))


if __name__ == "__main__":
    main()
