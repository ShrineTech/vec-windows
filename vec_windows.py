#!/usr/bin/env python3
"""vec-windows: Virtual Embryo Challenge Rules s.10 ("protected windows") as code.

Given the per-cell stage (and optionally genotype) labels of an EXTERNAL dataset you want to train on,
report which cells must be removed for a given task and why, citing the rule text, and write a
disclosure snippet for your method summary.

Stdlib only for CSV/TSV obs tables; reading .h5ad needs `anndata` (obs only, backed mode).

    python vec_windows.py check obs.csv --stage-col stage --task T1
    python vec_windows.py check atlas.h5ad --stage-col stage --task T2-heart --task T2-embryo \
        --genotype-col genotype --drop-out drop.csv --disclosure disclosure.md

Exit code: 0 = nothing to drop, 4 = some cells must be dropped (report printed), 1 = usage error.

THIS IS NOT AN OFFICIAL TOOL. The Rules are the authority, and they say: "If a source sits near the
boundary and you are unsure, ask before submitting." Rule text quoted from
https://virtualembryo.ai/challenge/rules s.10, read 2026-10-19.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, OrderedDict
from dataclasses import dataclass
from pathlib import Path

__version__ = "0.1.0"
EPS = 1e-9

# ---------------------------------------------------------------------------------------------------
# Rule text (verbatim, Rules s.10, read 2026-10-19)
# ---------------------------------------------------------------------------------------------------
RULE = {
    "extrap": ("External data is excluded from after E9.5 up to and including E13.5. E10.5 and E12.5, the "
               "validation and test targets themselves, are excluded absolutely. Anything after E13.5 may be "
               "used, provided the source and its stages are stated explicitly with the submission."),
    "extrap_ends": ("External data measured at exactly E9.5 is permitted [...] External data measured at "
                    "exactly E13.5 is not permitted"),
    "interp": ("external data may not be used if it stages strictly between the two released stages that "
               "bracket a held-out one. In Task 2 that means heart data strictly between E8.25 and E8.75, "
               "and embryo data strictly between E7.25 and E8.0."),
    "fractional": ("Fractional stagings do arise in external data, from somite counts and Theiler stages, and "
                   "they are placed on this scale and judged by the same boundary"),
    "genotype": ("Data from the held-out knockouts, from another allele of the same gene at a comparable "
                 "stage, or from a perturbation that phenocopies them, is treated as the held-out condition."),
    "general": ("A general-purpose resource that happens to span a held-out stage is not banned outright: the "
                "cells inside the window must be removed before training, and the removal described in the "
                "method summary."),
    "disclose": ("Every external source must be disclosed with the submission. An undisclosed source is a "
                 "rules violation whether or not it changed the result"),
    "ask": "If a source sits near the boundary and you are unsure, ask before submitting.",
}

# ---------------------------------------------------------------------------------------------------
# Task windows. Each window is (lo, hi, lo_inclusive, hi_inclusive, rule_key, label).
# ---------------------------------------------------------------------------------------------------
TASKS = {
    "T1": {"released": [8.5, 9.5], "held_out": [10.5, 12.5],
           "windows": [(9.5, 13.5, False, True, "extrap", "T1 extrapolation window (E9.5, E13.5]")],
           "disclose_after": 13.5},
    "T2-heart": {"released": [8.25, 8.75, 9.5], "held_out": [8.5, 10.5, 12.5],
                 "windows": [(8.25, 8.75, False, False, "interp", "T2 heart interpolation window (E8.25, E8.75)"),
                             (9.5, 13.5, False, True, "extrap", "T2 heart extrapolation window (E9.5, E13.5]")],
                 "disclose_after": 13.5},
    "T2-embryo": {"released": [6.75, 7.25, 8.0], "held_out": [7.5, 7.75],
                  "windows": [(7.25, 8.0, False, False, "interp", "T2 embryo interpolation window (E7.25, E8.0)")],
                  "disclose_after": None},
    "T3": {"released": [], "held_out": [], "windows": [], "disclose_after": None},
}
# Held-out genotypes (T3): Gata4 and beta-catenin (Ctnnb1) knockouts at E8.75. Matched case-insensitively.
HELD_OUT_GENOTYPES = {"Gata4": [r"gata[\s_-]?4"],
                      "beta-catenin (Ctnnb1)": [r"ctnnb1", r"(beta|β)[\s_-]?catenin", r"b[\s_-]?catenin"]}
WT_PATTERNS = [r"^\s*(wt|wild[\s_-]?type|control|ctrl|het|\+/\+)\s*$"]

CROSS_TASK_NOTE = {
    "T1": "Data-page vs Rules: none for T1 (external data is expressly encouraged).",
    "other": ("OPEN QUESTION: the Data page says 'Task 1 is the only task where training on external public "
              "single-cell data is allowed', while Rules s.10 permits external data for every task. Ask the "
              "organisers before using external data for T2/T3."),
}

# ---------------------------------------------------------------------------------------------------
# Stage parsing. Somite counts and Theiler stages map to an E-day RANGE (conservative). APPROXIMATE.
# Nominal dpc and somite ranges follow the standard Theiler (1989) staging as tabulated by eMouseAtlas
# (https://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/stagecriteria.html).
# UNVERIFIED in this release (the site was unreachable when this was written): every Theiler stage is
# widened to nominal +/- 0.5 day, so a label is only KEPT when its whole band is clear of a window.
# Override with --stage-map if you have a better staging for your source.
# ---------------------------------------------------------------------------------------------------
THEILER_NOMINAL_DPC = {10: 7.0, 11: 7.5, 12: 8.0, 13: 8.5, 14: 9.0, 15: 9.5, 16: 10.0, 17: 10.5, 18: 11.0,
                       19: 11.5, 20: 12.0, 21: 12.5, 22: 13.5, 23: 15.0, 24: 16.0, 25: 17.0, 26: 18.0}
THEILER_SOMITES = {12: (1, 7), 13: (8, 12), 14: (13, 20), 15: (21, 29), 16: (30, 34), 17: (35, 39),
                   18: (40, 44), 19: (45, 47), 20: (48, 51), 21: (52, 55)}
THEILER_BAND = 0.5

_E_RE = re.compile(r"^\s*E?\s*(\d+(?:\.\d+)?)\s*(?:dpc|d)?\s*$", re.I)
_TS_RE = re.compile(r"^\s*(?:TS|Theiler(?:\s*stage)?)\s*(\d+)\s*$", re.I)
_SOM_RE = re.compile(r"^\s*(?:(\d+)\s*(?:s|ss|som|somites?|-?somite)|(?:ss|somites?)\s*(\d+))\s*$", re.I)


@dataclass(frozen=True)
class Stage:
    lo: float
    hi: float
    how: str          # "E-day", "Theiler", "somite", "map", "unknown"

    @property
    def point(self):
        return abs(self.hi - self.lo) < EPS


def theiler_range(ts: int) -> tuple:
    d = THEILER_NOMINAL_DPC[ts]
    return d - THEILER_BAND, d + THEILER_BAND


def parse_stage(label, stage_map: dict | None = None) -> Stage:
    s = "" if label is None else str(label).strip()
    if stage_map and s in stage_map:
        lo, hi = stage_map[s]
        return Stage(lo, hi, "map")
    m = _E_RE.match(s)
    if m:
        v = float(m.group(1))
        return Stage(v, v, "E-day")
    m = _TS_RE.match(s)
    if m and int(m.group(1)) in THEILER_NOMINAL_DPC:
        return Stage(*theiler_range(int(m.group(1))), "Theiler")
    m = _SOM_RE.match(s)
    if m:
        n = int(m.group(1) or m.group(2))
        for ts, (a, b) in THEILER_SOMITES.items():
            if a <= n <= b:
                return Stage(*theiler_range(ts), "somite")
    return Stage(float("nan"), float("nan"), "unknown")


def _in_window(x: float, w) -> bool:
    lo, hi, lo_inc, hi_inc = w[:4]
    above = x > lo - EPS if lo_inc else x > lo + EPS
    below = x < hi + EPS if hi_inc else x < hi - EPS
    return above and below


def _range_hits(st: Stage, w) -> str:
    """'all' / 'some' / 'none' of the stage range [lo, hi] falls inside window w."""
    if st.point:
        return "all" if _in_window(st.lo, w) else "none"
    lo, hi = w[0], w[1]
    probes = [st.lo, st.hi] + [v for v in (lo, hi) if st.lo - EPS <= v <= st.hi + EPS] + \
             [v + d for v in (lo, hi) for d in (1e-6, -1e-6) if st.lo <= v + d <= st.hi]
    hits = [_in_window(p, w) for p in probes]
    return "all" if all(hits) else ("some" if any(hits) else "none")


def genotype_hit(g) -> str | None:
    if g is None:
        return None
    s = str(g)
    if any(re.search(p, s, re.I) for p in WT_PATTERNS):
        return None
    for name, pats in HELD_OUT_GENOTYPES.items():
        if any(re.search(p, s, re.I) for p in pats):
            return name
    return None


def verdict(label, tasks: list, genotype=None, stage_map=None) -> tuple:
    """-> (action, reason, rule_key). action in KEEP, KEEP_DISCLOSE, DROP."""
    st = parse_stage(label, stage_map)
    g = genotype_hit(genotype)
    if g:
        return ("DROP", f"genotype {genotype!r} matches held-out knockout {g} (Task 3; dropped for every task: "
                        "'comparable stage' is not defined numerically, so this is conservative)", "genotype")
    if st.how == "unknown":
        return ("DROP", f"stage {label!r} could not be placed on the E-day scale; dropped until staged "
                        "(use --stage-map)", "fractional")
    for t in tasks:
        for w in TASKS[t]["windows"]:
            h = _range_hits(st, w)
            where = f"E{st.lo:g}" if st.point else f"E{st.lo:g}-E{st.hi:g} ({st.how}, approximate)"
            if h == "all":
                absolute = any(abs(st.lo - x) < EPS for x in TASKS[t]["held_out"]) and st.point
                tag = " -- a held-out target itself, excluded absolutely" if absolute else ""
                return ("DROP", f"{where} is inside the {w[5]}{tag}", w[4])
            if h == "some":
                return ("DROP", f"{where} straddles the {w[5]} boundary; dropped conservatively. "
                                + RULE["ask"], w[4])
    for t in tasks:
        da = TASKS[t]["disclose_after"]
        if da is not None and st.lo > da + EPS:
            return ("KEEP_DISCLOSE", f"after E{da:g}: allowed for {t}, but the source and its stages must be "
                                     "stated explicitly with the submission", "extrap")
    return ("KEEP", "outside every protected window for " + ", ".join(tasks), "disclose")


# ---------------------------------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------------------------------
def read_obs(path: Path, cols: list) -> list:
    """Return a list of dicts with only `cols` (missing columns -> error)."""
    if path.suffix == ".h5ad":
        try:
            import anndata as ad
        except ImportError:
            raise SystemExit("reading .h5ad needs anndata (pip install anndata), or export obs to CSV first")
        a = ad.read_h5ad(path, backed="r")
        obs = a.obs
        miss = [c for c in cols if c not in obs.columns]
        if miss:
            raise SystemExit(f"{path}: obs has no column(s) {miss}; have {list(obs.columns)[:20]}")
        out = [{"cell": str(i), **{c: obs[c].iloc[k] for c in cols}} for k, i in enumerate(obs.index)]
        a.file.close()
        return out
    delim = "\t" if path.suffix in (".tsv", ".tab") else ","
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter=delim)
        miss = [c for c in cols if c not in (r.fieldnames or [])]
        if miss:
            raise SystemExit(f"{path}: no column(s) {miss}; have {r.fieldnames}")
        first = (r.fieldnames or [""])[0]
        return [{"cell": row.get("cell") or row.get(first), **{c: row[c] for c in cols}} for row in r]


def read_stage_map(path) -> dict:
    """CSV with columns label,lo,hi (E-days)."""
    with open(path, newline="") as fh:
        return {r["label"].strip(): (float(r["lo"]), float(r["hi"])) for r in csv.DictReader(fh)}


def check(rows: list, stage_col: str, tasks: list, genotype_col=None, stage_map=None) -> dict:
    per_label: "OrderedDict[tuple, dict]" = OrderedDict()
    actions = Counter()
    drop_cells = []
    for r in rows:
        lab, gt = r[stage_col], (r.get(genotype_col) if genotype_col else None)
        key = (str(lab), None if gt is None else str(gt))
        if key not in per_label:
            a, why, rk = verdict(lab, tasks, gt, stage_map)
            per_label[key] = {"stage": key[0], "genotype": key[1], "action": a, "reason": why, "rule": rk, "n": 0}
        e = per_label[key]
        e["n"] += 1
        actions[e["action"]] += 1
        if e["action"] == "DROP":
            drop_cells.append((r["cell"], e["stage"], e["reason"]))
    return {"tasks": tasks, "n_cells": len(rows), "actions": dict(actions), "labels": list(per_label.values()),
            "drop_cells": drop_cells}


def disclosure(res: dict, source: str, stage_col: str) -> str:
    kept = res["n_cells"] - res["actions"].get("DROP", 0)
    lab = lambda l: l["stage"] + (f" [{l['genotype']}]" if l["genotype"] not in (None, "") else "")
    dropped_labels = sorted({lab(l) for l in res["labels"] if l["action"] == "DROP"})
    kept_labels = sorted({lab(l) for l in res["labels"] if l["action"] != "DROP"})
    late = sorted({l["stage"] for l in res["labels"] if l["action"] == "KEEP_DISCLOSE"})
    lines = [
        f"External data: {source}.",
        f"Stages used (column `{stage_col}`): {', '.join(kept_labels) or 'none'}.",
        f"Removed before training, per Rules s.10 protected windows for {', '.join(res['tasks'])}: "
        f"{res['actions'].get('DROP', 0)} of {res['n_cells']} cells, stage/genotype labels "
        f"{', '.join(dropped_labels) or 'none'}. {kept} cells retained.",
    ]
    if late:
        lines.append(f"Stages after E13.5 used, stated explicitly as required: {', '.join(late)}.")
    lines.append("Filtering done with vec-windows " + __version__ + " (somite/Theiler stagings mapped to "
                 "conservative approximate E-day ranges; boundary cases dropped).")
    return "\n".join(lines)


def format_report(res: dict) -> str:
    out = [f"vec-windows {__version__}  tasks={','.join(res['tasks'])}  cells={res['n_cells']}  "
           f"actions={res['actions']}"]
    w = max([len(l["stage"]) for l in res["labels"]] + [5])
    for l in res["labels"]:
        gt = f" [{l['genotype']}]" if l["genotype"] not in (None, "") else ""
        out.append(f"  {l['action']:<13} {l['stage']:<{w}}{gt}  n={l['n']:<7} {l['reason']}")
    used = sorted({l["rule"] for l in res["labels"] if l["action"] != "KEEP"})
    if used:
        out.append("\nRule text (Rules s.10, verbatim):")
        out += [f"  [{k}] \"{RULE[k]}\"" for k in used]
    out.append("\n  [general] \"" + RULE["general"] + "\"")
    out.append("  [disclose] \"" + RULE["disclose"] + "\"")
    if any(t != "T1" for t in res["tasks"]):
        out.append("\n" + CROSS_TASK_NOTE["other"])
    out.append("NOT OFFICIAL. Genotype phenocopies cannot be detected from labels; review your source. "
               + RULE["ask"])
    return "\n".join(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="vec-windows", description="Rules s.10 protected-window checker")
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("obs", help=".csv/.tsv obs table (first column = cell id) or .h5ad")
    c.add_argument("--stage-col", required=True)
    c.add_argument("--task", action="append", required=True, choices=list(TASKS))
    c.add_argument("--genotype-col")
    c.add_argument("--stage-map", help="CSV label,lo,hi to override stage parsing")
    c.add_argument("--source", default="<name, version, URL, licence>")
    c.add_argument("--drop-out", help="write dropped cell ids + reason as CSV")
    c.add_argument("--disclosure", help="write a method-summary disclosure snippet (markdown)")
    a = p.parse_args(argv)
    cols = [a.stage_col] + ([a.genotype_col] if a.genotype_col else [])
    rows = read_obs(Path(a.obs), cols)
    res = check(rows, a.stage_col, a.task, a.genotype_col, read_stage_map(a.stage_map) if a.stage_map else None)
    print(format_report(res))
    if a.drop_out:
        with open(a.drop_out, "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["cell", "stage", "reason"])
            wr.writerows(res["drop_cells"])
    snip = disclosure(res, a.source, a.stage_col)
    if a.disclosure:
        Path(a.disclosure).write_text(snip + "\n")
    print("\nDisclosure snippet:\n" + snip)
    return 4 if res["actions"].get("DROP") else 0


if __name__ == "__main__":
    sys.exit(main())
