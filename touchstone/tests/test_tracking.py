from collections import Counter

from touchstone.tracking import candidate_metrics, stack_rows


def _result(name, consensus, reward, stack, stress=None):
    r = {
        "structure": f"/designs/{name}.pdb",
        "coordination_number": 4,
        "donors": ["N", "N", "N", "S"],
        "consensus": consensus,
        "reward": reward,
        "stack": stack,
    }
    if stress is not None:
        r["stress"] = stress
    return r


_RAN = lambda stage, label, score: {"stage": stage, "status": "ran", "label": label, "score": score}  # noqa: E731


class TestStackRows:
    def test_columns_and_cells_track_the_stack_view(self):
        stack = [_RAN("geometry", "trust", 0.27), _RAN("mlip", "weak", 0.565),
                 {"stage": "mogul", "status": "needs_input", "detail": "a CSD licence"}]
        cols, rows, counts = stack_rows([_result("ni_motif_02", "weak", 0.41, stack)])

        assert cols == ["design", "CN", "donors", "geometry", "mlip", "mogul", "consensus", "reward"]
        row = rows[0]
        assert row[0] == "ni_motif_02.pdb" and row[2] == "NNNS"
        assert row[3] == "trust 0.27"  # ran ⇒ "label score"
        assert row[4] == "weak 0.56"
        assert row[5] == "needs_input"  # not run ⇒ status
        assert row[-2:] == ["weak", 0.41]  # consensus + reward passthrough
        assert counts == Counter({"weak": 1})

    def test_stress_holds_column_counts_non_defer_conditions(self):
        stack = [_RAN("geometry", "trust", 0.27)]
        stress = {"neutral": {"label": "trust"}, "leachate": {"label": "weak"}, "low_pH": {"label": "defer"}}
        cols, rows, _ = stack_rows([_result("d", "weak", 0.4, stack, stress)])
        assert "stress_holds" in cols
        assert rows[0][cols.index("stress_holds")] == "2/3"  # trust + weak hold; low_pH defers

    def test_no_stress_column_when_absent(self):
        cols, _, counts = stack_rows([
            _result("a", "trust", 0.9, [_RAN("geometry", "trust", 0.9)]),
            _result("b", "defer", 0.0, [_RAN("geometry", "defer", 0.1)]),
        ])
        assert "stress_holds" not in cols
        assert counts == Counter({"trust": 1, "defer": 1})


class TestCandidateMetrics:
    def test_per_model_scores_and_consensus_weight(self):
        stack = [_RAN("geometry", "trust", 0.27), _RAN("mlip", "weak", 0.57),
                 {"stage": "mogul", "status": "needs_input", "detail": "..."}]
        m = candidate_metrics(_result("d", "weak", 0.41, stack))
        assert m["reward"] == 0.41
        assert m["consensus_weight"] == 0.5  # weak → 0.5
        assert m["score/geometry"] == 0.27 and m["score/mlip"] == 0.57
        assert "score/mogul" not in m  # only verifiers that actually ran get a series
