"""Tests for touchstone.search — novelty screening over sequence and structure corpora.

The subprocess calls are not exercised; the parsing and ranking are, because that is where
both known bugs lived. Ranking sequence hits by maximum percent identity once reported 88.8%
for a design whose genuine best hit was 36% over 83% coverage, and a FoldSeek line that fails
to parse must drop rather than take the whole batch with it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from touchstone.search import (
    MIN_QUERY_COVERAGE,
    FoldHit,
    SequenceHit,
    SequenceSearch,
    kmer_corpus,
    parse_foldseek,
)


def tabular(*rows: tuple[str, str, float, int, int, float, float, float]) -> str:
    """An MMseqs2 --format-output block: query,target,pident,alnlen,qlen,qcov,evalue,bits."""
    return "\n".join("\t".join(str(field) for field in row) for row in rows)


def search(**kwargs) -> SequenceSearch:
    return SequenceSearch(binary=Path("mmseqs"), databases=(), **kwargs)


class TestSequenceRanking:
    """Bitscore above a coverage floor, never maximum identity."""

    def test_prefers_the_higher_bitscore_not_the_higher_identity(self):
        """The 88.8%-vs-36% bug: a short high-identity alignment must not win."""
        best: dict[str, SequenceHit] = {}
        search().absorb(
            tabular(
                (
                    "design",
                    "sliver",
                    88.8,
                    9,
                    90,
                    0.9,
                    1e-2,
                    20.0,
                ),  # high identity, low bits
                ("design", "real_hit", 36.0, 75, 90, 0.83, 1e-20, 150.0),
            ),
            best,
        )
        assert best["design"].target == "real_hit"
        assert best["design"].identity_pct == 36.0

    def test_drops_hits_below_the_coverage_floor(self):
        best: dict[str, SequenceHit] = {}
        search().absorb(
            tabular(("design", "tiny", 99.0, 5, 90, 0.1, 1e-3, 500.0)), best
        )
        assert best == {}

    def test_coverage_floor_is_configurable(self):
        best: dict[str, SequenceHit] = {}
        search(min_coverage=0.05).absorb(
            tabular(("design", "tiny", 99.0, 5, 90, 0.1, 1e-3, 500.0)), best
        )
        assert best["design"].target == "tiny"

    def test_identity_given_as_a_fraction_is_scaled(self):
        """MMseqs reports pident as a fraction or a percentage depending on version."""
        best: dict[str, SequenceHit] = {}
        search().absorb(
            tabular(("design", "hit", 0.36, 75, 90, 0.83, 1e-20, 150.0)), best
        )
        assert best["design"].identity_pct == pytest.approx(36.0)

    def test_a_percentage_is_left_alone(self):
        best: dict[str, SequenceHit] = {}
        search().absorb(
            tabular(("design", "hit", 36.0, 75, 90, 0.83, 1e-20, 150.0)), best
        )
        assert best["design"].identity_pct == pytest.approx(36.0)

    def test_merges_across_corpora_keeping_the_best(self):
        """Two databases are searched in turn into the same dict."""
        best: dict[str, SequenceHit] = {}
        engine = search()
        engine.absorb(
            tabular(("d", "swissprot_hit", 40.0, 70, 90, 0.8, 1e-9, 100.0)), best
        )
        engine.absorb(tabular(("d", "pdb_hit", 38.0, 72, 90, 0.8, 1e-30, 200.0)), best)
        assert best["d"].target == "pdb_hit"

    def test_separates_queries(self):
        best: dict[str, SequenceHit] = {}
        search().absorb(
            tabular(
                ("a", "ha", 40.0, 70, 90, 0.8, 1e-9, 100.0),
                ("b", "hb", 50.0, 70, 90, 0.8, 1e-9, 120.0),
            ),
            best,
        )
        assert set(best) == {"a", "b"}

    @pytest.mark.parametrize(
        "line", ["", "too\tfew\tfields", "d\tt\tnot_a_number\t1\t1\t0.9\t1\t1"]
    )
    def test_malformed_lines_are_skipped_not_raised(self, line):
        best: dict[str, SequenceHit] = {}
        search().absorb(line, best)
        assert best == {}

    def test_the_coverage_floor_default_is_half(self):
        assert MIN_QUERY_COVERAGE == 0.5


class TestParseFoldseek:
    """FoldSeek tabular output -> hits per design."""

    def test_reads_tm_and_coverage(self):
        hits = parse_foldseek(
            "design_binder.pdb\t2znvF00\t0.893\t0.888\t0.9\t1e-10\t200"
        )
        assert hits["design"] == [FoldHit(0.893, 0.888, "2znvF00")]

    @pytest.mark.parametrize(
        "query,stem",
        [
            ("d_binder.pdb", "d"),
            ("d.pdb", "d"),
            ("d.cif", "d"),
            ("d", "d"),
        ],
    )
    def test_strips_query_suffixes(self, query, stem):
        assert stem in parse_foldseek(f"{query}\tt\t0.5\t0.9\t0.9\t1e-5\t100")

    def test_collects_multiple_hits_per_design(self):
        hits = parse_foldseek(
            "\n".join(
                [
                    "d.pdb\tfirst\t0.9\t0.9\t0.9\t1e-10\t200",
                    "d.pdb\tsecond\t0.6\t0.8\t0.8\t1e-5\t100",
                ]
            )
        )
        assert [h.target for h in hits["d"]] == ["first", "second"]

    @pytest.mark.parametrize(
        "text", ["", "short\tline", "d\tt\tNaNscore\t0.9\t0.9\t1\t1"]
    )
    def test_malformed_output_yields_nothing_rather_than_raising(self, text):
        assert parse_foldseek(text) == {} or all(
            v == [] for v in parse_foldseek(text).values()
        )

    def test_one_bad_line_does_not_lose_the_good_ones(self):
        hits = parse_foldseek(
            "\n".join(
                [
                    "d.pdb\tgood\t0.9\t0.9\t0.9\t1e-10\t200",
                    "garbage",
                ]
            )
        )
        assert [h.target for h in hits["d"]] == ["good"]


class TestKmerCorpus:
    """Exact k-mer screening, decisive where an alignment score is arguable."""

    def test_finds_an_exact_kmer_present_in_the_corpus(self, tmp_path):
        flat = tmp_path / "corpus.txt"
        flat.write_text("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGR\n")
        assert kmer_corpus([flat])(["MQIFVKTLTGKTITLEVEP"]) == {"MQIFVKTLTGKTITLEVEP"}

    def test_reports_nothing_for_a_sequence_not_present(self, tmp_path):
        flat = tmp_path / "corpus.txt"
        flat.write_text("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ\n")
        assert kmer_corpus([flat])(["WWWWWWWWWWWWWWWWWWWW"]) == set()

    def test_no_kmers_means_no_corpus_scan(self, tmp_path):
        assert kmer_corpus([tmp_path / "absent.txt"])([]) == set()

    def test_searches_every_corpus(self, tmp_path):
        first, second = tmp_path / "a.txt", tmp_path / "b.txt"
        first.write_text("AAAAAAAAAA\n")
        second.write_text("CCCCCCCCCC\n")
        assert kmer_corpus([first, second])(["AAAAAAAAAA", "CCCCCCCCCC"]) == {
            "AAAAAAAAAA",
            "CCCCCCCCCC",
        }
