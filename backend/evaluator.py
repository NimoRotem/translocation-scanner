#!/usr/bin/env python3
"""Post-run evaluator — standalone, never called by the pipeline.

Usage:
    python evaluator.py results/calls.json ground_truth.json [--window 1000000]

Loads ground truth events and matches pipeline calls against them.
Reports: TP, FP, FN, sensitivity, precision, rank of each truth event.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass


@dataclass
class TruthEvent:
    chrom_a: str
    pos_a: int
    chrom_b: str
    pos_b: int
    label: str = ""
    match_window: int = 1_000_000


@dataclass
class PipelineCall:
    chrom_a: str
    pos_a: int
    chrom_b: str
    pos_b: int
    tier: str = ""
    score: float = 0.0
    cluster_id: str = ""
    rank: int = 0


def load_ground_truth(path: str) -> list[TruthEvent]:
    with open(path) as f:
        data = json.load(f)
    events = []
    for item in data.get("events", []):
        events.append(TruthEvent(
            chrom_a=item["chrom_a"],
            pos_a=item["pos_a"],
            chrom_b=item["chrom_b"],
            pos_b=item["pos_b"],
            label=item.get("label", ""),
            match_window=item.get("match_window", 1_000_000),
        ))
    return events


def load_pipeline_calls(path: str) -> list[PipelineCall]:
    with open(path) as f:
        data = json.load(f)

    # Handle both direct list and nested structure
    if isinstance(data, list):
        calls_raw = data
    elif isinstance(data, dict):
        calls_raw = data.get("calls", data.get("validated_calls", []))
    else:
        calls_raw = []

    calls = []
    for rank, item in enumerate(calls_raw, 1):
        calls.append(PipelineCall(
            chrom_a=item.get("chrom_a", ""),
            pos_a=item.get("pos_a", 0),
            chrom_b=item.get("chrom_b", ""),
            pos_b=item.get("pos_b", 0),
            tier=item.get("tier", ""),
            score=item.get("score", 0),
            cluster_id=item.get("cluster_id", ""),
            rank=rank,
        ))
    return calls


def _normalize_chrom(c: str) -> str:
    """Strip 'chr' prefix for matching."""
    return c.replace("chr", "")


def match_events(
    truth: list[TruthEvent],
    calls: list[PipelineCall],
    default_window: int = 1_000_000,
) -> dict:
    """Match pipeline calls to truth events within window."""
    matched_truth: dict[int, PipelineCall | None] = {i: None for i in range(len(truth))}
    matched_calls: set[int] = set()

    # For each truth event, find the best matching call
    for ti, te in enumerate(truth):
        window = te.match_window or default_window
        best_call = None
        best_dist = float('inf')

        for ci, call in enumerate(calls):
            if ci in matched_calls:
                continue

            # Check both orientations of the match
            match_ab = (
                _normalize_chrom(call.chrom_a) == _normalize_chrom(te.chrom_a)
                and _normalize_chrom(call.chrom_b) == _normalize_chrom(te.chrom_b)
                and abs(call.pos_a - te.pos_a) <= window
                and abs(call.pos_b - te.pos_b) <= window
            )
            match_ba = (
                _normalize_chrom(call.chrom_a) == _normalize_chrom(te.chrom_b)
                and _normalize_chrom(call.chrom_b) == _normalize_chrom(te.chrom_a)
                and abs(call.pos_a - te.pos_b) <= window
                and abs(call.pos_b - te.pos_a) <= window
            )

            if match_ab or match_ba:
                dist = abs(call.pos_a - te.pos_a) + abs(call.pos_b - te.pos_b)
                if dist < best_dist:
                    best_dist = dist
                    best_call = ci

        if best_call is not None:
            matched_truth[ti] = calls[best_call]
            matched_calls.add(best_call)

    # Compute metrics
    tp = sum(1 for v in matched_truth.values() if v is not None)
    fn = sum(1 for v in matched_truth.values() if v is None)
    fp = len(calls) - len(matched_calls)

    sensitivity = tp / len(truth) if truth else 0
    precision = tp / len(calls) if calls else 0

    # Truth event details
    truth_details = []
    for ti, te in enumerate(truth):
        match = matched_truth[ti]
        truth_details.append({
            "truth_event": {
                "chrom_a": te.chrom_a, "pos_a": te.pos_a,
                "chrom_b": te.chrom_b, "pos_b": te.pos_b,
                "label": te.label,
            },
            "matched": match is not None,
            "matched_call": {
                "cluster_id": match.cluster_id,
                "rank": match.rank,
                "tier": match.tier,
                "score": match.score,
            } if match else None,
        })

    return {
        "summary": {
            "total_truth": len(truth),
            "total_calls": len(calls),
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "sensitivity": round(sensitivity, 4),
            "precision": round(precision, 4),
        },
        "truth_events": truth_details,
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate pipeline calls against ground truth")
    parser.add_argument("calls", help="Pipeline calls JSON file")
    parser.add_argument("ground_truth", help="Ground truth JSON file")
    parser.add_argument("--window", type=int, default=1_000_000,
                        help="Match window in bp (default: 1Mb)")
    parser.add_argument("--json", action="store_true",
                        help="Output raw JSON instead of summary")
    args = parser.parse_args()

    truth = load_ground_truth(args.ground_truth)
    calls = load_pipeline_calls(args.calls)

    results = match_events(truth, calls, default_window=args.window)

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        s = results["summary"]
        print(f"Ground truth events: {s['total_truth']}")
        print(f"Pipeline calls:      {s['total_calls']}")
        print(f"True positives:      {s['tp']}")
        print(f"False positives:     {s['fp']}")
        print(f"False negatives:     {s['fn']}")
        print(f"Sensitivity:         {s['sensitivity']:.1%}")
        print(f"Precision:           {s['precision']:.1%}")
        print()

        for td in results["truth_events"]:
            te = td["truth_event"]
            status = "FOUND" if td["matched"] else "MISSED"
            rank_info = ""
            if td["matched_call"]:
                mc = td["matched_call"]
                rank_info = f" -> rank {mc['rank']}, tier={mc['tier']}, score={mc['score']}"
            print(f"  [{status}] {te['chrom_a']}:{te['pos_a']} ↔ "
                  f"{te['chrom_b']}:{te['pos_b']} ({te['label']}){rank_info}")


if __name__ == "__main__":
    main()
