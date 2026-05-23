from __future__ import annotations

import ast
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _assert_close(name: str, actual: float, expected: float, rel_tol: float = 1.0e-9) -> None:
    if not math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=0.0):
        raise AssertionError(f"{name}: expected {expected}, got {actual}")


def _load_source_constants() -> tuple[dict[str, dict[str, float]], dict[str, float]]:
    source_path = ROOT / "ni604_opb_fatigue.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))

    studless_scf = None
    sn_defaults: dict[str, float] = {}

    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "STUDLESS_SCF":
                    studless_scf = ast.literal_eval(node.value)
        elif isinstance(node, ast.ClassDef) and node.name == "SNCurve":
            for item in node.body:
                if isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
                    if item.target.id in {"log10_k", "m"} and item.value is not None:
                        sn_defaults[item.target.id] = float(ast.literal_eval(item.value))

    if studless_scf is None:
        raise AssertionError("Could not find STUDLESS_SCF in ni604_opb_fatigue.py")
    if {"log10_k", "m"} - set(sn_defaults):
        raise AssertionError("Could not find SNCurve default S-N parameters")

    return studless_scf, sn_defaults


def check_code_constants() -> None:
    studless_scf, sn_defaults = _load_source_constants()
    expected_scf = {
        "A": {"TT": 4.48, "OPB": 0.00, "IPB": 1.25},
        "B": {"TT": 2.08, "OPB": 1.06, "IPB": 0.71},
        "C": {"TT": 1.04, "OPB": 1.21, "IPB": 1.50},
    }

    for hotspot, modes in expected_scf.items():
        for mode, expected in modes.items():
            _assert_close(f"SCF {hotspot} {mode}", studless_scf[hotspot][mode], expected)

    _assert_close("Free-corrosion S-N log10(K)", sn_defaults["log10_k"], 12.436)
    _assert_close("Free-corrosion S-N slope m", sn_defaults["m"], 3.0)


def check_table8_summary() -> None:
    rows = _read_csv(ROOT / "validation" / "mdpi_2024_table8_damage_life.csv")
    by_quantity = {row["quantity"]: row for row in rows}

    damage = by_quantity["fatigue_damage"]
    damage_values = [float(damage[name]) for name in ["TT", "Hotspot_A", "Hotspot_B", "Hotspot_C"]]
    if not damage_values[0] < damage_values[1] < damage_values[2] < damage_values[3]:
        raise AssertionError("Table 8 damage ordering should be TT < A < B < C")

    life = by_quantity["fatigue_life_years"]
    life_values = [float(life[name]) for name in ["TT", "Hotspot_A", "Hotspot_B", "Hotspot_C"]]
    if not life_values[0] > life_values[1] > life_values[2] > life_values[3]:
        raise AssertionError("Table 8 life ordering should be TT > A > B > C")

    safety_factor = 10.0
    for name in ["TT", "Hotspot_A", "Hotspot_B", "Hotspot_C"]:
        expected_life = 1.0 / (float(damage[name]) * safety_factor)
        published_life = float(life[name])
        _assert_close(
            f"Table 8 life relation {name}",
            published_life,
            expected_life,
            rel_tol=6.0e-3,
        )


def check_appendix_a4_subset() -> None:
    rows = _read_csv(ROOT / "validation" / "mdpi_2024_appendix_a4_subset.csv")
    hotspot_columns = ["Hotspot_A", "Hotspot_B", "Hotspot_C"]

    governing_counts: dict[str, int] = {}
    for row in rows:
        values = {name: float(row[name]) for name in hotspot_columns}
        actual = max(values, key=values.get)
        expected = row["expected_governing_hotspot"]
        if actual != expected:
            raise AssertionError(f"Case {row['case']}: expected {expected}, got {actual}")
        governing_counts[actual] = governing_counts.get(actual, 0) + 1

    if governing_counts.get("Hotspot_A", 0) == 0 or governing_counts.get("Hotspot_C", 0) == 0:
        raise AssertionError("Benchmark subset should include both tension- and OPB-governed cases")


def main() -> None:
    check_code_constants()
    check_table8_summary()
    check_appendix_a4_subset()
    print("MDPI 2024 literature benchmark checks passed.")


if __name__ == "__main__":
    main()
