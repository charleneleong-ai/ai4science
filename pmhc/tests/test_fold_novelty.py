"""The fold gate's verdict, which a TM-score alone cannot give.

Every one of the first 11 designs matched something at TM 0.71-0.93, so the numbers below are
the real ones from that batch: a 95.1%-identity match to 5xs3 chain B is a known protein handed
back, the 0.0%-identity match to 7zwj chain A at a higher TM is not, and the immunoglobulin
beta-sandwich hit is the class the organisers' own check has already failed us on.
"""

from __future__ import annotations

import pytest
from fold_novelty import FoldHit, best_hits


def hit(target: str, tm: float, identity: float, design: str = "d") -> FoldHit:
    return FoldHit(
        design=design,
        database="pdb",
        target=target,
        identity=identity,
        tm_score=tm,
        evalue=1e-9,
        length=80,
    )


class TestVerdict:
    @pytest.mark.parametrize(
        "target,tm,identity,expected",
        [
            ("5xs3-assembly1_B", 1.005, 0.951, "copy_of_known"),
            ("7zwj-assembly1_A", 0.995, 0.000, "shared_fold"),
            ("af_O97394_1370_1467_2.60.40.10", 0.715, 0.026, "ig_sandwich"),
            ("af_O94651_10_125_1.20.58.70", 0.870, 0.035, "shared_fold"),
            ("1hfeT00", 0.410, 0.020, "distinct_fold"),
        ],
    )
    def test_reads_tm_and_identity_together(self, target, tm, identity, expected):
        assert hit(target, tm, identity).verdict == expected

    def test_a_copied_immunoglobulin_reads_as_a_copy(self):
        """copy_of_known outranks ig_sandwich, since it is the stronger statement."""
        assert hit("af_X_2.60.40.10", 0.9, 0.95).verdict == "copy_of_known"


class TestBestHits:
    def write(self, tmp_path, rows):
        table = tmp_path / "hits.m8"
        table.write_text(
            "".join(
                f"{d}\t{t}\t{i}\t{tm}\t1e-9\t80\n" for d, t, i, tm in rows
            )
        )
        return table

    def test_the_strongest_match_defines_the_fold(self, tmp_path):
        """Not the most disqualifying one. A permissive hit list holds a 2.60.40 entry for
        almost every design, so reading its tail reports noise: ranking that way called 8 of
        11 designs immunoglobulin, one of them over its own TM 0.865 hit to another fold."""
        table = self.write(
            tmp_path,
            [
                ("design_a", "af_Q_1.20.58.70", 0.03, 0.93),
                ("design_a", "af_O97394_1370_1467_2.60.40.10", 0.03, 0.54),
            ],
        )
        [found] = best_hits(table, "cath50").values()
        assert found.verdict == "shared_fold"
        assert found.tm_score == 0.93

    def test_an_immunoglobulin_top_hit_is_still_flagged(self, tmp_path):
        table = self.write(
            tmp_path,
            [
                ("design_a", "af_O97394_1370_1467_2.60.40.10", 0.03, 0.72),
                ("design_a", "af_Q_1.20.58.70", 0.03, 0.44),
            ],
        )
        [found] = best_hits(table, "cath50").values()
        assert found.verdict == "ig_sandwich"

    def test_skips_short_lines(self, tmp_path):
        table = tmp_path / "hits.m8"
        table.write_text("design_a\taf_Q_1.20.58.70\n\ndesign_a\tt\t0.0\t0.9\t1e-9\t80\n")
        assert len(best_hits(table, "pdb")) == 1
