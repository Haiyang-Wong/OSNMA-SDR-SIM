#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict


def subframe_start_tow(tow, start_time):
    page_num = ((tow - start_time) // 2) % 15
    return tow - page_num * 2


def parse_args():
    parser = argparse.ArgumentParser(
        description="Count generated cross-auth tags per subframe from osnma_tag_log.csv."
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        default="osnma_tag_log.csv",
        help="Path to osnma_tag_log.csv. Defaults to ./osnma_tag_log.csv.",
    )
    parser.add_argument(
        "--start-time",
        type=int,
        default=303061,
        help="Simulation start_time used to derive the 15-page subframe start TOW.",
    )
    parser.add_argument(
        "--time-field",
        choices=("write_tow", "auth_tow"),
        default="write_tow",
        help="TOW field used for grouping. Defaults to write_tow.",
    )
    parser.add_argument(
        "--nonzero-only",
        action="store_true",
        help="Print only subframes with at least one cross tag.",
    )
    parser.add_argument(
        "--show-pairs",
        action="store_true",
        help="Also print PRNa->PRNd pair counts for each subframe.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    counts = Counter()
    pairs_by_subframe = defaultdict(Counter)
    rows_read = 0
    cross_rows = 0

    with open(args.csv_path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        required = {"write_tow", "auth_tow", "prna", "prnd"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Missing required CSV columns: {', '.join(sorted(missing))}")

        for row in reader:
            rows_read += 1
            prna = int(row["prna"])
            prnd = int(row["prnd"])
            if prna == prnd:
                continue

            cross_rows += 1
            tow = int(row[args.time_field])
            subframe_tow = subframe_start_tow(tow, args.start_time)
            counts[subframe_tow] += 1
            pairs_by_subframe[subframe_tow][(prna, prnd)] += 1

    print(f"CSV: {args.csv_path}")
    print(f"time_field: {args.time_field}")
    print(f"rows_read: {rows_read}")
    print(f"cross_tags_prna_ne_prnd: {cross_rows}")
    print()
    if args.show_pairs:
        print("subframe_tow,cross_tag_count,pairs")
    else:
        print("subframe_tow,cross_tag_count")

    for subframe_tow in sorted(counts):
        count = counts[subframe_tow]
        if args.nonzero_only and count == 0:
            continue
        if args.show_pairs:
            pair_text = " ".join(
                f"{prna}->{prnd}:{pair_count}"
                for (prna, prnd), pair_count in pairs_by_subframe[subframe_tow].most_common()
            )
            print(f"{subframe_tow},{count},{pair_text}")
        else:
            print(f"{subframe_tow},{count}")


if __name__ == "__main__":
    main()
