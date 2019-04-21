import sys


def main():
    with open(sys.argv[1]) as f:
        for line in f:
            line = line.strip().split()
            if line[0] == "S":
                print(line[3])


if __name__ == "__main__":
    main()
