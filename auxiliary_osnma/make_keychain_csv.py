#!/usr/bin/env python3
"""Build a compact OSNMA keychain CSV from full I/NAV page dumps.

Output format:

    TOW,OSNMA40

The extractor first scans all raw pages and measures OSNMA40 coverage per PRN.
It then builds one keychain stream by keeping the current source PRN while it
has valid OSNMA40 at the current TOW. When that source disappears, the next
source is the PRN with the longest continuous OSNMA40 coverage from that TOW.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_INPUT = Path(__file__).with_name("13_MAY_2026_GST_12_11_01.csv")
DEFAULT_OUTPUT = Path(__file__).with_name("13_MAY_2026_GST_12_11_01_keychain.csv")
ZERO_OSNMA40 = "0" * 40


@dataclass
class PrnCoverage:
    prn: int
    valid_rows: int = 0
    min_tow: int | None = None
    max_tow: int | None = None
    zero_osnma_rows: int = 0
    duplicate_tows: int = 0
    first_seen_order: int = 0
    tows: set[int] = field(default_factory=set)
    segments: list[tuple[int, int, int, int]] = field(default_factory=list)


@dataclass
class ScanResult:
    tow_prn_osnma: dict[int, dict[int, str]]
    prn_coverage: dict[int, PrnCoverage]
    input_rows: int
    valid_rows: int
    invalid_rows: int
    empty_rows: int
    malformed_rows: int
    hex_parse_error_rows: int
    short_page_rows: int
    duplicate_tow_prn_rows: int
    min_tow: int | None
    max_tow: int | None
    global_tows: list[int]


def hex_to_bits(value: str) -> str:
    value = value.strip()

    if value.lower().startswith("0x"):
        value = value[2:]

    if len(value) % 2:
        raise ValueError("hex string has odd length")

    raw = bytes.fromhex(value)
    return "".join(f"{byte:08b}" for byte in raw)


def report_path_for(output_path: Path) -> Path:
    return output_path.with_name(f"{output_path.stem}_report.csv")


def split_segments(tows: list[int], expected_step: int) -> list[tuple[int, int, int, int]]:
    if not tows:
        return []

    segments: list[tuple[int, int, int, int]] = []
    start = tows[0]
    previous = tows[0]
    count = 1

    for tow in tows[1:]:
        if tow - previous == expected_step:
            count += 1
        else:
            segments.append((start, previous, previous - start, count))
            start = tow
            count = 1
        previous = tow

    segments.append((start, previous, previous - start, count))
    return segments


def scan_input(input_path: Path, expected_step: int) -> ScanResult:
    tow_prn_osnma: dict[int, dict[int, str]] = defaultdict(dict)
    prn_coverage: dict[int, PrnCoverage] = {}
    prn_first_seen: dict[int, int] = {}

    input_rows = 0
    valid_rows = 0
    invalid_rows = 0
    empty_rows = 0
    malformed_rows = 0
    hex_parse_error_rows = 0
    short_page_rows = 0
    duplicate_tow_prn_rows = 0

    with input_path.open(newline="") as src:
        reader = csv.reader(src)

        for row in reader:
            input_rows += 1

            if not row or all(not col.strip() for col in row):
                empty_rows += 1
                invalid_rows += 1
                continue

            if len(row) < 4:
                malformed_rows += 1
                invalid_rows += 1
                continue

            try:
                tow = int(row[0].strip())
                prn = int(row[2].strip())
            except (ValueError, TypeError):
                malformed_rows += 1
                invalid_rows += 1
                continue

            try:
                page_bits = hex_to_bits(row[3].strip())
            except (ValueError, TypeError):
                hex_parse_error_rows += 1
                invalid_rows += 1
                continue

            if len(page_bits) < 178:
                short_page_rows += 1
                invalid_rows += 1
                continue

            osnma40 = page_bits[138:178]
            if len(osnma40) != 40:
                short_page_rows += 1
                invalid_rows += 1
                continue

            if prn not in prn_first_seen:
                prn_first_seen[prn] = len(prn_first_seen)

            coverage = prn_coverage.setdefault(
                prn,
                PrnCoverage(prn=prn, first_seen_order=prn_first_seen[prn]),
            )
            coverage.valid_rows += 1
            coverage.min_tow = tow if coverage.min_tow is None else min(coverage.min_tow, tow)
            coverage.max_tow = tow if coverage.max_tow is None else max(coverage.max_tow, tow)

            if osnma40 == ZERO_OSNMA40:
                coverage.zero_osnma_rows += 1

            if prn in tow_prn_osnma[tow]:
                duplicate_tow_prn_rows += 1
                coverage.duplicate_tows += 1
                continue

            tow_prn_osnma[tow][prn] = osnma40
            coverage.tows.add(tow)
            valid_rows += 1

    for coverage in prn_coverage.values():
        coverage.segments = split_segments(sorted(coverage.tows), expected_step)

    global_tows = sorted(tow_prn_osnma)

    return ScanResult(
        tow_prn_osnma=dict(tow_prn_osnma),
        prn_coverage=prn_coverage,
        input_rows=input_rows,
        valid_rows=valid_rows,
        invalid_rows=invalid_rows,
        empty_rows=empty_rows,
        malformed_rows=malformed_rows,
        hex_parse_error_rows=hex_parse_error_rows,
        short_page_rows=short_page_rows,
        duplicate_tow_prn_rows=duplicate_tow_prn_rows,
        min_tow=global_tows[0] if global_tows else None,
        max_tow=global_tows[-1] if global_tows else None,
        global_tows=global_tows,
    )


def coverage_length(
    prn: int,
    current_tow: int,
    tow_prn_osnma: dict[int, dict[int, str]],
    expected_step: int,
) -> int:
    length = 0
    tow = current_tow

    while tow_prn_osnma.get(tow, {}).get(prn, ZERO_OSNMA40) != ZERO_OSNMA40:
        length += 1
        tow += expected_step

    return length


def select_source_prn(
    current_tow: int,
    tow_prn_osnma: dict[int, dict[int, str]],
    prn_coverage: dict[int, PrnCoverage],
    expected_step: int,
) -> tuple[int | None, int]:
    candidates = tow_prn_osnma.get(current_tow, {})
    if not candidates:
        return None, 0

    ranked: list[tuple[int, int, int, int]] = []
    for prn, osnma40 in candidates.items():
        if osnma40 == ZERO_OSNMA40:
            continue
        length = coverage_length(prn, current_tow, tow_prn_osnma, expected_step)
        first_seen_order = prn_coverage[prn].first_seen_order
        ranked.append((-length, prn, first_seen_order, prn))

    if not ranked:
        return None, 0

    best = min(ranked)
    return best[3], -best[0]


def next_source_tow(
    source_prn: int,
    current_tow: int,
    tow_prn_osnma: dict[int, dict[int, str]],
    prn_coverage: dict[int, PrnCoverage],
) -> int | None:
    for tow in sorted(prn_coverage[source_prn].tows):
        if tow >= current_tow and tow_prn_osnma.get(tow, {}).get(source_prn) != ZERO_OSNMA40:
            return tow
    return None


def build_expected_tows(global_tows: list[int], expected_step: int) -> list[int]:
    if not global_tows:
        return []

    return list(range(global_tows[0], global_tows[-1] + expected_step, expected_step))


def write_outputs(
    output_path: Path,
    report_path: Path,
    rows: list[tuple[int, str, int | None, bool, str]],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as dst:
        writer = csv.writer(dst, lineterminator="\n")
        for tow, osnma40, _source_prn, _is_zero_fill, _reason in rows:
            writer.writerow([tow, osnma40])

    with report_path.open("w", newline="") as dst:
        writer = csv.writer(dst, lineterminator="\n")
        writer.writerow(["TOW", "OSNMA40", "source_prn", "is_zero_fill", "reason"])
        for tow, osnma40, source_prn, is_zero_fill, reason in rows:
            writer.writerow(
                [
                    tow,
                    osnma40,
                    "" if source_prn is None else source_prn,
                    str(is_zero_fill).lower(),
                    reason,
                ]
            )


def build_keychain(
    input_path: Path,
    output_path: Path,
    expected_step: int,
    max_source_gap_fill: int,
) -> dict[str, object]:
    scan = scan_input(input_path, expected_step)
    expected_tows = build_expected_tows(scan.global_tows, expected_step)

    current_source_prn: int | None = None
    rows: list[tuple[int, str, int | None, bool, str]] = []
    source_switch_count = 0
    zero_fill_tows = 0
    source_usage: Counter[int] = Counter()
    source_minmax: dict[int, list[int]] = {}
    source_segments_raw: list[tuple[int, int, int]] = []

    current_segment_prn: int | None = None
    current_segment_start: int | None = None
    current_segment_end: int | None = None

    def close_segment() -> None:
        nonlocal current_segment_prn, current_segment_start, current_segment_end
        if current_segment_prn is not None and current_segment_start is not None:
            source_segments_raw.append(
                (current_segment_prn, current_segment_start, current_segment_end or current_segment_start)
            )
        current_segment_prn = None
        current_segment_start = None
        current_segment_end = None

    for tow in expected_tows:
        available = scan.tow_prn_osnma.get(tow, {})

        current_source_osnma = available.get(current_source_prn) if current_source_prn is not None else None

        if current_source_prn is not None and current_source_osnma not in (None, ZERO_OSNMA40):
            osnma40 = available[current_source_prn]
            reason = "use_existing_source"
            source_prn = current_source_prn
            is_zero_fill = False
        elif current_source_prn is not None:
            next_tow = next_source_tow(
                current_source_prn,
                tow,
                scan.tow_prn_osnma,
                scan.prn_coverage,
            )
            if next_tow is not None and next_tow - tow <= max_source_gap_fill:
                source_prn, _span = select_source_prn(
                    tow,
                    scan.tow_prn_osnma,
                    scan.prn_coverage,
                    expected_step,
                )

                if source_prn is None:
                    osnma40 = ZERO_OSNMA40
                    reason = "source_gap_zero_fill"
                    source_prn = current_source_prn
                    is_zero_fill = True
                    zero_fill_tows += 1
                else:
                    osnma40 = available[source_prn]
                    reason = "temporary_source_fill"
                    is_zero_fill = False
            else:
                previous_source_prn = current_source_prn
                source_prn, span = select_source_prn(
                    tow,
                    scan.tow_prn_osnma,
                    scan.prn_coverage,
                    expected_step,
                )

                if source_prn is None:
                    osnma40 = ZERO_OSNMA40
                    reason = "no_available_osnma"
                    is_zero_fill = True
                    zero_fill_tows += 1
                    current_source_prn = None
                    close_segment()
                else:
                    osnma40 = available[source_prn]
                    reason = "switch_source"
                    is_zero_fill = False
                    current_source_prn = source_prn
                    source_switch_count += 1
                    print(
                        "source_switch: "
                        f"tow={tow}, "
                        f"from={previous_source_prn}, "
                        f"to={source_prn}, "
                        f"continuous_rows={span}",
                        flush=True,
                    )
        else:
            previous_source_prn = current_source_prn
            source_prn, span = select_source_prn(
                tow,
                scan.tow_prn_osnma,
                scan.prn_coverage,
                expected_step,
            )

            if source_prn is None:
                osnma40 = ZERO_OSNMA40
                reason = "no_available_osnma"
                is_zero_fill = True
                zero_fill_tows += 1
                current_source_prn = None
                close_segment()
            else:
                osnma40 = available[source_prn]
                reason = "switch_source"
                is_zero_fill = False
                current_source_prn = source_prn
                source_switch_count += 1
                print(
                    "source_switch: "
                    f"tow={tow}, "
                    f"from={'' if previous_source_prn is None else previous_source_prn}, "
                    f"to={source_prn}, "
                    f"continuous_rows={span}",
                    flush=True,
                )

        rows.append((tow, osnma40, source_prn, is_zero_fill, reason))

        if source_prn is None:
            continue

        source_usage[source_prn] += 1
        if source_prn not in source_minmax:
            source_minmax[source_prn] = [tow, tow]
        else:
            source_minmax[source_prn][1] = tow

        if current_segment_prn != source_prn:
            close_segment()
            current_segment_prn = source_prn
            current_segment_start = tow
        current_segment_end = tow

    close_segment()

    report_path = report_path_for(output_path)
    write_outputs(output_path, report_path, rows)

    return {
        "scan": scan,
        "input_rows": scan.input_rows,
        "valid_rows": scan.valid_rows,
        "invalid_rows": scan.invalid_rows,
        "empty_rows": scan.empty_rows,
        "malformed_rows": scan.malformed_rows,
        "hex_parse_error_rows": scan.hex_parse_error_rows,
        "short_page_rows": scan.short_page_rows,
        "duplicate_tow_prn_rows": scan.duplicate_tow_prn_rows,
        "output_tows": len(rows),
        "min_tow": expected_tows[0] if expected_tows else None,
        "max_tow": expected_tows[-1] if expected_tows else None,
        "zero_fill_tows": zero_fill_tows,
        "source_switch_count": source_switch_count,
        "source_usage": source_usage,
        "source_minmax": source_minmax,
        "source_prn_segments": source_segments_raw,
        "output_path": str(output_path),
        "report_path": str(report_path),
    }


def read_output_tows(output_path: Path) -> list[int]:
    tows: list[int] = []

    if not output_path.exists():
        return tows

    with output_path.open(newline="") as src:
        reader = csv.reader(src)

        for row in reader:
            if not row:
                continue

            try:
                tows.append(int(row[0].strip()))
            except (ValueError, TypeError):
                continue

    return tows


def check_tow_continuity(
    output_path: Path,
    expected_step: int = 2,
) -> dict[str, int | list[tuple[int, int, int, list[int]]]]:
    tows = sorted(read_output_tows(output_path))
    gaps: list[tuple[int, int, int, list[int]]] = []

    for current_tow, next_tow in zip(tows, tows[1:]):
        diff = next_tow - current_tow

        if diff != expected_step:
            missing_tows = list(range(current_tow + expected_step, next_tow, expected_step))
            gaps.append((current_tow, next_tow, diff, missing_tows))

    return {
        "checked_tows": len(tows),
        "expected_step": expected_step,
        "missing_segments": len(gaps),
        "gaps": gaps,
    }


def print_scan_summary(scan: ScanResult) -> None:
    print("prn_coverage_summary:")
    for prn in sorted(scan.prn_coverage):
        coverage = scan.prn_coverage[prn]
        print(
            f"  PRN{prn}: "
            f"valid_osnma_rows={coverage.valid_rows}, "
            f"unique_tows={len(coverage.tows)}, "
            f"min_tow={coverage.min_tow}, "
            f"max_tow={coverage.max_tow}, "
            f"segments={len(coverage.segments)}, "
            f"zero_osnma_rows={coverage.zero_osnma_rows}, "
            f"has_zero_osnma={coverage.zero_osnma_rows > 0}, "
            f"duplicate_tows={coverage.duplicate_tows}"
        )
        for index, (start, end, duration, count) in enumerate(coverage.segments[:8], start=1):
            print(
                f"    segment{index}: "
                f"start_tow={start}, end_tow={end}, duration={duration}, count={count}"
            )
        hidden = len(coverage.segments) - 8
        if hidden > 0:
            print(f"    additional_segments_not_shown={hidden}")


def print_stats(stats: dict[str, object]) -> None:
    for key in (
        "input_rows",
        "valid_rows",
        "invalid_rows",
        "empty_rows",
        "malformed_rows",
        "hex_parse_error_rows",
        "short_page_rows",
        "duplicate_tow_prn_rows",
        "output_tows",
        "min_tow",
        "max_tow",
        "zero_fill_tows",
        "source_switch_count",
        "output_path",
        "report_path",
    ):
        print(f"{key}={stats[key]}")

    print("source_prn_segments:")
    for prn, start, end in stats["source_prn_segments"]:
        print(f"  PRN{prn}: start_tow={start}, end_tow={end}, duration={end - start}")

    print("source_prn_usage:")
    source_usage: Counter[int] = stats["source_usage"]
    source_minmax: dict[int, list[int]] = stats["source_minmax"]
    for prn, count in sorted(source_usage.items()):
        start, end = source_minmax[prn]
        print(f"  PRN{prn}: used_tows={count}, start_tow={start}, end_tow={end}")

    print_scan_summary(stats["scan"])


def print_output_status(output_path: Path) -> None:
    exists = output_path.exists()
    size = output_path.stat().st_size if exists else 0

    print(f"output={output_path}")
    print(f"output_exists={exists}")
    print(f"output_size={size}")


def print_tow_continuity(output_path: Path, expected_step: int) -> None:
    continuity = check_tow_continuity(output_path, expected_step)
    gaps = continuity["gaps"]
    checked_tows = continuity["checked_tows"]

    if checked_tows == 0:
        print("tow_continuity=NO_DATA")
        print(f"expected_step={expected_step}")
        print("missing_segments=0")
        return

    if not gaps:
        print("tow_continuity=OK")
        print(f"tow_step={expected_step}")
        print("missing_segments=0")
        return

    first_gap = gaps[0]

    print("tow_continuity=FAILED")
    print(f"expected_step={expected_step}")
    print(f"missing_segments={continuity['missing_segments']}")
    print(f"first_missing_after_tow={first_gap[0]}")
    print(f"first_bad_next_tow={first_gap[1]}")

    for current_tow, next_tow, diff, missing_tows in gaps[:20]:
        print(
            "gap: "
            f"current_tow={current_tow}, "
            f"next_tow={next_tow}, "
            f"diff={diff}, "
            f"missing_tows={missing_tows}"
        )

    hidden = len(gaps) - 20
    if hidden > 0:
        print(f"additional_gaps_not_shown={hidden}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input raw CSV. Format: TOW,GST_WN,PRN,full_nav_page_hex",
    )

    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output compact CSV. Format: TOW,OSNMA40",
    )

    parser.add_argument(
        "--expected-step",
        type=int,
        default=2,
        help="Expected TOW step in seconds. Default: 2",
    )

    parser.add_argument(
        "--max-source-gap-fill",
        type=int,
        default=30,
        help=(
            "Keep the current source PRN across gaps up to this many seconds by "
            "writing zero OSNMA40. Default: 30"
        ),
    )

    args = parser.parse_args()

    print(f"input_path={args.input}", flush=True)
    print(f"output_path={args.output}", flush=True)
    print("source_prn_policy=dynamic_longest_continuous_coverage", flush=True)

    if not args.input.exists():
        print(f"error=input file does not exist: {args.input}", file=sys.stderr)
        return 2

    stats = build_keychain(
        input_path=args.input,
        output_path=args.output,
        expected_step=args.expected_step,
        max_source_gap_fill=args.max_source_gap_fill,
    )

    print_stats(stats)
    print_output_status(args.output)
    print_tow_continuity(args.output, args.expected_step)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
