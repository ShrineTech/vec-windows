"""Tests for vec_windows on SYNTHETIC label tables only. Stdlib unittest; pytest also collects it.

    python -m unittest -v test_vec_windows
"""
import csv
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import vec_windows as vw  # noqa: E402


def act(label, *tasks, genotype=None):
    return vw.verdict(label, list(tasks), genotype)[0]


class TestT1Extrapolation(unittest.TestCase):
    def test_boundaries(self):
        cases = {"E8.5": "KEEP", "E9.5": "KEEP", "9.5": "KEEP", "E9.51": "DROP", "E9.75": "DROP",
                 "E10.5": "DROP", "E12.5": "DROP", "E13.5": "DROP", "E13.6": "KEEP_DISCLOSE",
                 "E14.5": "KEEP_DISCLOSE", "E7.75": "KEEP", "E6.5": "KEEP"}
        for lab, want in cases.items():
            self.assertEqual(act(lab, "T1"), want, lab)

    def test_target_marked_absolute(self):
        self.assertIn("excluded absolutely", vw.verdict("E10.5", ["T1"])[1])
        self.assertNotIn("excluded absolutely", vw.verdict("E11.0", ["T1"])[1])


class TestT2Interpolation(unittest.TestCase):
    def test_heart(self):
        cases = {"E8.25": "KEEP", "E8.3": "DROP", "E8.5": "DROP", "E8.74": "DROP", "E8.75": "KEEP",
                 "E9.0": "KEEP", "E9.5": "KEEP", "E10.0": "DROP", "E13.5": "DROP", "E14.0": "KEEP_DISCLOSE"}
        for lab, want in cases.items():
            self.assertEqual(act(lab, "T2-heart"), want, lab)

    def test_embryo(self):
        cases = {"E6.75": "KEEP", "E7.25": "KEEP", "E7.3": "DROP", "E7.5": "DROP", "E7.75": "DROP",
                 "E7.99": "DROP", "E8.0": "KEEP", "E8.5": "KEEP", "E10.5": "KEEP"}
        for lab, want in cases.items():
            self.assertEqual(act(lab, "T2-embryo"), want, lab)

    def test_union_of_tasks(self):
        self.assertEqual(act("E7.5", "T1"), "KEEP")
        self.assertEqual(act("E7.5", "T1", "T2-embryo"), "DROP")
        self.assertEqual(act("E8.5", "T1", "T2-heart"), "DROP")


class TestFractionalStagings(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(vw.parse_stage("TS17").how, "Theiler")
        self.assertEqual(vw.parse_stage("Theiler stage 15").how, "Theiler")
        self.assertEqual(vw.parse_stage("25 somites").how, "somite")
        self.assertEqual(vw.parse_stage("ss25").how, "somite")
        self.assertEqual(vw.parse_stage("E9.5").lo, 9.5)
        self.assertEqual(vw.parse_stage("adult").how, "unknown")

    def test_theiler_and_somites_conservative(self):
        self.assertEqual(act("TS13", "T1"), "KEEP")          # ~E8.5 +/- 0.5, clear of (9.5, 13.5]
        self.assertEqual(act("TS15", "T1"), "DROP")          # ~E9.5 band straddles the boundary
        self.assertIn("straddles", vw.verdict("TS15", ["T1"])[1])
        self.assertEqual(act("TS17", "T1"), "DROP")          # ~E10.5 fully inside
        self.assertEqual(act("TS23", "T1"), "KEEP_DISCLOSE")  # ~E15
        self.assertEqual(act("10 somites", "T2-heart"), "DROP")   # TS13 ~E8.5 band hits (8.25, 8.75)
        self.assertEqual(act("4 somites", "T2-embryo"), "DROP")   # TS12 ~E8.0 band hits (7.25, 8.0)

    def test_unknown_stage_dropped(self):
        self.assertEqual(act("P0", "T1"), "DROP")
        self.assertEqual(act("", "T1"), "DROP")

    def test_stage_map_override(self):
        self.assertEqual(vw.verdict("batch_A", ["T1"], None, {"batch_A": (8.9, 9.2)})[0], "KEEP")


class TestGenotype(unittest.TestCase):
    def test_held_out_knockouts(self):
        for g in ["Gata4 KO", "gata4-/-", "Gata4_cKO", "Ctnnb1 flox/flox", "beta-catenin cKO", "β-catenin KO"]:
            self.assertEqual(act("E8.75", "T3", genotype=g), "DROP", g)
        for g in ["WT", "wild-type", "Mab21l2 KO", "control", None]:
            self.assertEqual(act("E8.75", "T3", genotype=g), "KEEP", g)


class TestEndToEnd(unittest.TestCase):
    def _write_synthetic(self, d):
        rows = [("cell_%d" % i, s, g) for i, (s, g) in enumerate(
            [("E9.5", "WT")] * 5 + [("E10.5", "WT")] * 3 + [("E14.5", "WT")] * 2 +
            [("TS15", "WT")] * 4 + [("E8.75", "Gata4 KO")] * 1)]
        p = Path(d) / "obs.csv"
        with p.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cell", "stage", "genotype"])
            w.writerows(rows)
        return p

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write_synthetic(d)
            out, drop, disc = io.StringIO(), Path(d) / "drop.csv", Path(d) / "disc.md"
            with redirect_stdout(out):
                rc = vw.main(["check", str(p), "--stage-col", "stage", "--genotype-col", "genotype",
                              "--task", "T1", "--drop-out", str(drop), "--disclosure", str(disc),
                              "--source", "SYNTHETIC test atlas v0"])
            self.assertEqual(rc, 4)
            dropped = list(csv.DictReader(drop.open()))
            self.assertEqual(len(dropped), 3 + 4 + 1)
            self.assertEqual({r["stage"] for r in dropped}, {"E10.5", "TS15", "E8.75"})
            text = disc.read_text()
            self.assertIn("8 of 15 cells", text)
            self.assertIn("E14.5", text)                       # after-E13.5 stated explicitly
            self.assertIn("excluded from after E9.5 up to and including E13.5", out.getvalue())

    def test_h5ad_if_anndata(self):
        try:
            import anndata as ad
            import numpy as np
            import pandas as pd
        except ImportError:
            self.skipTest("anndata not installed")
        with tempfile.TemporaryDirectory() as d:
            obs = pd.DataFrame({"stage": ["E9.5", "E11.5", "E15.5"]}, index=["a", "b", "c"])
            p = Path(d) / "syn.h5ad"
            ad.AnnData(X=np.zeros((3, 2), np.float32), obs=obs).write_h5ad(p)
            res = vw.check(vw.read_obs(p, ["stage"]), "stage", ["T1"])
            self.assertEqual(res["actions"], {"KEEP": 1, "DROP": 1, "KEEP_DISCLOSE": 1})
            self.assertEqual(res["drop_cells"][0][0], "b")


if __name__ == "__main__":
    unittest.main()
