"""
BV NI604-style top-chain OPB/IPB combined fatigue helper.

The calculation follows the public NI604 workflow:
    tension / OPB moment / IPB moment time series
    -> hotspot stress time series
    -> rainflow counted stress ranges
    -> S-N damage with Miner summation

Important limitation:
    OPB fatigue cannot be derived from mooring tension alone. For OPB/IPB
    combined fatigue, provide either interlink moments or interlink angles
    from a top-chain / fairlead local model. If only tension is available,
    this module calculates the TT contribution with zero bending moments.

Units used by this module:
    tension: kN
    chain diameter: mm
    moment: N.mm
    stress: MPa = N/mm^2
"""

from __future__ import annotations

import argparse
import json
import math
from collections import deque
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


@dataclass(frozen=True)
class ChainConfig:
    """Top-chain data for the fatigue calculation."""

    diameter_mm: float
    pretension_kN: float
    mbl_kN: float | None = None
    chain_type: str = "studless"
    z_corr: float = 1.08
    z_stiffness: float = 1.06
    friction_coefficient: float = 0.30
    gamma_tt_override: float | None = None

    @property
    def ipb_alpha(self) -> float:
        chain_type = self.chain_type.lower()
        if chain_type == "studless":
            return 2.33
        if chain_type == "studlink":
            return 2.06
        raise ValueError("chain_type must be 'studless' or 'studlink'")

    @property
    def gamma_tt(self) -> float:
        if self.gamma_tt_override is not None:
            return self.gamma_tt_override
        if self.mbl_kN is None or self.mbl_kN <= 0.0:
            raise ValueError("mbl_kN or gamma_tt_override is required for hotspot C")
        return max(1.0 + 0.9 * (self.pretension_kN / self.mbl_kN - 0.15), 0.95)


@dataclass(frozen=True)
class SNCurve:
    """S-N curve in the form N = 10**log10_k / stress_range**m."""

    log10_k: float = 12.436
    m: float = 3.0

    def cycles_to_failure(self, stress_range_MPa: np.ndarray) -> np.ndarray:
        ranges = np.asarray(stress_range_MPa, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.power(10.0, self.log10_k) / np.power(ranges, self.m)


STUDLESS_SCF = {
    "A": {"TT": 4.48, "OPB": 0.00, "IPB": 1.25},
    "B": {"TT": 2.08, "OPB": 1.06, "IPB": 0.71},
    "B_PRIME": {"TT": 1.65, "OPB": 1.15, "IPB": 0.66},
    "C": {"TT": 1.04, "OPB": 1.21, "IPB": 1.50},
}


def ni604_interlink_moment_Nmm(
    interlink_angle_deg: Iterable[float] | np.ndarray,
    tension_kN: Iterable[float] | np.ndarray,
    diameter_mm: float,
    friction_coefficient: float = 0.30,
    apply_sliding_limit: bool = True,
) -> np.ndarray:
    """Evaluate the NI604 Appendix 1 parametric interlink moment model.

    The formula is intended for chain diameters from 84 mm to 146 mm and
    chain grades covered by NI604 Appendix 1. Outside that range, treat
    results as a screening value requiring class/project confirmation.
    """

    angle = np.asarray(interlink_angle_deg, dtype=float)
    tension = np.asarray(tension_kN, dtype=float)
    abs_angle = np.abs(angle)

    # Appendix 1 expresses the nonlinear interlink bending moment with
    # empirical angle functions P(alpha), a(alpha), and b(alpha). The input
    # angle is kept signed, but the empirical formula uses its magnitude.
    c = 354.0
    g = 0.93
    p_alpha = abs_angle + 0.307 * abs_angle**3 + 0.048 * abs_angle**5
    a_alpha = 0.439 + 0.532 * np.tanh(1.020 * abs_angle)
    b_alpha = 0.433 + 1.640 * np.tanh(1.320 * abs_angle)

    normalized_tension = tension / (0.14 * diameter_mm**2)
    diameter_factor = diameter_mm / 100.0

    with np.errstate(divide="ignore", invalid="ignore"):
        moment_abs = (
            math.pi
            * diameter_mm**3
            / 16.0
            * c
            * p_alpha
            / (g + p_alpha)
            * np.power(np.maximum(normalized_tension, 0.0), a_alpha)
            * np.power(diameter_factor, 2.0 * a_alpha + b_alpha)
        )

    moment = np.sign(angle) * np.nan_to_num(
        moment_abs,
        nan=0.0,
        posinf=np.finfo(float).max,
        neginf=0.0,
    )

    if apply_sliding_limit:
        # Sliding at the contact limits the transferable interlink moment.
        # Keep the sign from the prescribed interlink rotation.
        threshold = friction_coefficient * tension * 1000.0 * diameter_mm / 2.0
        moment = np.sign(moment) * np.minimum(np.abs(moment), np.maximum(threshold, 0.0))

    return moment


def _as_array(values: Iterable[float] | np.ndarray, n: int | None = None) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if n is not None and arr.size != n:
        raise ValueError(f"Expected {n} values, got {arr.size}")
    return arr


def stress_components_MPa(
    tension_kN: Iterable[float] | np.ndarray,
    opb_moment_Nmm: Iterable[float] | np.ndarray,
    ipb_moment_Nmm: Iterable[float] | np.ndarray,
    config: ChainConfig,
    hotspot: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return TT, OPB and IPB hotspot stress time series in MPa."""

    hotspot_key = hotspot.upper()
    if hotspot_key not in STUDLESS_SCF:
        raise ValueError(f"Unknown hotspot '{hotspot}'. Use one of {sorted(STUDLESS_SCF)}")

    tension = _as_array(tension_kN)
    opb_moment = _as_array(opb_moment_Nmm, tension.size)
    ipb_moment = _as_array(ipb_moment_Nmm, tension.size)
    d = config.diameter_mm

    scf = dict(STUDLESS_SCF[hotspot_key])
    if hotspot_key == "C":
        scf["OPB"] *= config.gamma_tt

    # Stress components are computed as time series so rainflow can be applied
    # after TT, OPB, and IPB are combined at each hotspot.
    sigma_tt = scf["TT"] * 2.0 * tension * 1000.0 / (math.pi * d**2)
    sigma_opb = scf["OPB"] * 16.0 * opb_moment / (math.pi * d**3)
    sigma_ipb = scf["IPB"] * config.ipb_alpha * ipb_moment / (math.pi * d**3)
    return sigma_tt, sigma_opb, sigma_ipb


def combined_stress_MPa(
    tension_kN: Iterable[float] | np.ndarray,
    opb_moment_Nmm: Iterable[float] | np.ndarray,
    ipb_moment_Nmm: Iterable[float] | np.ndarray,
    config: ChainConfig,
    hotspot: str,
    opb_sign: int = 1,
    ipb_sign: int = 1,
) -> np.ndarray:
    """Return combined hotspot stress time series for one OPB/IPB sign choice."""

    sigma_tt, sigma_opb, sigma_ipb = stress_components_MPa(
        tension_kN=tension_kN,
        opb_moment_Nmm=opb_moment_Nmm,
        ipb_moment_Nmm=ipb_moment_Nmm,
        config=config,
        hotspot=hotspot,
    )
    # The OPB/IPB signs depend on the physical hotspot side. For screening,
    # fatigue_from_dataframe evaluates all sign combinations and governs by
    # maximum damage.
    return config.z_corr * (
        sigma_tt
        + opb_sign * config.z_stiffness * sigma_opb
        + ipb_sign * config.z_stiffness * sigma_ipb
    )


def _reversals(series: Iterable[float] | np.ndarray) -> list[tuple[int, float]]:
    values = np.asarray(series, dtype=float)
    if values.size == 0:
        return []

    mask = np.isfinite(values)
    if not np.all(mask):
        values = values[mask]
    if values.size <= 1:
        return [(0, float(values[0]))] if values.size else []

    compressed: list[tuple[int, float]] = []
    for i, value in enumerate(values):
        if not compressed or value != compressed[-1][1]:
            compressed.append((i, float(value)))

    if len(compressed) <= 2:
        return compressed

    result = [compressed[0]]
    for previous, current, following in zip(compressed, compressed[1:], compressed[2:]):
        left = current[1] - previous[1]
        right = following[1] - current[1]
        if left * right < 0.0:
            result.append(current)
    result.append(compressed[-1])
    return result


def rainflow_cycles(series: Iterable[float] | np.ndarray) -> pd.DataFrame:
    """ASTM E1049-style rainflow counting.

    Returns one row per counted cycle with columns:
        range, mean, count, start_index, end_index
    """

    points: deque[tuple[int, float]] = deque()
    cycles: list[dict[str, float]] = []

    def add_cycle(p1: tuple[int, float], p2: tuple[int, float], count: float) -> None:
        stress_range = abs(p2[1] - p1[1])
        if stress_range <= 0.0:
            return
        cycles.append(
            {
                "range": stress_range,
                "mean": 0.5 * (p1[1] + p2[1]),
                "count": count,
                "start_index": p1[0],
                "end_index": p2[0],
            }
        )

    for point in _reversals(series):
        points.append(point)

        # ASTM-style stack counting: when the newest range closes or exceeds
        # the previous range, the older range becomes a counted cycle.
        while len(points) >= 3:
            x1, x2, x3 = points[-3][1], points[-2][1], points[-1][1]
            newer_range = abs(x3 - x2)
            older_range = abs(x2 - x1)

            if newer_range < older_range:
                break

            if len(points) == 3:
                add_cycle(points[0], points[1], 0.5)
                points.popleft()
            else:
                add_cycle(points[-3], points[-2], 1.0)
                last = points.pop()
                points.pop()
                points.pop()
                points.append(last)

    while len(points) > 1:
        add_cycle(points[0], points[1], 0.5)
        points.popleft()

    return pd.DataFrame(cycles, columns=["range", "mean", "count", "start_index", "end_index"])


def fatigue_damage_from_stress(
    stress_MPa: Iterable[float] | np.ndarray,
    sn_curve: SNCurve = SNCurve(),
    exposure_scale: float = 1.0,
) -> tuple[float, pd.DataFrame]:
    """Calculate Miner damage from a stress time series."""

    cycles = rainflow_cycles(stress_MPa)
    if cycles.empty:
        cycles["cycles_to_failure"] = []
        cycles["damage"] = []
        return 0.0, cycles

    ranges = cycles["range"].to_numpy(dtype=float)
    cycles_to_failure = sn_curve.cycles_to_failure(ranges)
    damage = cycles["count"].to_numpy(dtype=float) / cycles_to_failure
    cycles["cycles_to_failure"] = cycles_to_failure
    cycles["damage"] = damage * exposure_scale
    return float(cycles["damage"].sum()), cycles


def _duration_seconds(
    row_count: int,
    dt_s: float | None = None,
    time_s: Iterable[float] | np.ndarray | None = None,
) -> float | None:
    if time_s is not None:
        time = np.asarray(time_s, dtype=float)
        if time.size >= 2:
            return float(np.nanmax(time) - np.nanmin(time))
    if dt_s is not None and row_count >= 2:
        return float((row_count - 1) * dt_s)
    return None


def fatigue_from_dataframe(
    df: pd.DataFrame,
    config: ChainConfig,
    sn_curve: SNCurve = SNCurve(),
    tension_column: str = "tension_kN",
    opb_moment_column: str = "opb_moment_Nmm",
    ipb_moment_column: str = "ipb_moment_Nmm",
    opb_angle_column: str = "opb_interlink_angle_deg",
    ipb_angle_column: str = "ipb_interlink_angle_deg",
    time_column: str | None = None,
    dt_s: float | None = None,
    exposure_seconds: float | None = None,
    occurrence_probability: float = 1.0,
    hotspots: tuple[str, ...] = ("A", "B", "B_PRIME", "C"),
) -> dict[str, object]:
    """Calculate fatigue damage for each hotspot from an input dataframe.

    Required input column:
        tension_kN

    Optional bending input:
        opb_moment_Nmm / ipb_moment_Nmm
        or opb_interlink_angle_deg / ipb_interlink_angle_deg

    If moments are absent but angles are present, the NI604 Appendix 1
    moment model is used. If neither is present, zero bending is assumed.
    """

    if tension_column not in df.columns:
        raise ValueError(f"Missing required column: {tension_column}")

    tension = df[tension_column].to_numpy(dtype=float)
    n = tension.size

    if opb_moment_column in df.columns:
        opb_moment = df[opb_moment_column].to_numpy(dtype=float)
        opb_source = opb_moment_column
    elif opb_angle_column in df.columns:
        # If only interlink angles are available, derive moments from the
        # same nonlinear relation used by the FEM rotational springs.
        opb_moment = ni604_interlink_moment_Nmm(
            df[opb_angle_column].to_numpy(dtype=float),
            tension,
            config.diameter_mm,
            friction_coefficient=config.friction_coefficient,
        )
        opb_source = opb_angle_column
    else:
        opb_moment = np.zeros(n)
        opb_source = "zero"

    if ipb_moment_column in df.columns:
        ipb_moment = df[ipb_moment_column].to_numpy(dtype=float)
        ipb_source = ipb_moment_column
    elif ipb_angle_column in df.columns:
        ipb_moment = ni604_interlink_moment_Nmm(
            df[ipb_angle_column].to_numpy(dtype=float),
            tension,
            config.diameter_mm,
            friction_coefficient=config.friction_coefficient,
        )
        ipb_source = ipb_angle_column
    else:
        ipb_moment = np.zeros(n)
        ipb_source = "zero"

    duration = _duration_seconds(
        row_count=n,
        dt_s=dt_s,
        time_s=df[time_column].to_numpy(dtype=float) if time_column else None,
    )
    exposure_scale = occurrence_probability
    if exposure_seconds is not None:
        if duration is None or duration <= 0.0:
            raise ValueError("dt_s or time_column is required when exposure_seconds is used")
        # Scale rainflow damage from the simulated duration to the requested
        # exposure duration; occurrence_probability handles sea-state weighting.
        exposure_scale *= exposure_seconds / duration

    results: dict[str, object] = {
        "input": {
            "samples": int(n),
            "duration_seconds": duration,
            "exposure_scale": exposure_scale,
            "opb_source": opb_source,
            "ipb_source": ipb_source,
            "gamma_tt": config.gamma_tt,
        },
        "hotspots": {},
    }

    for hotspot in hotspots:
        sign_results = []
        for opb_sign in (-1, 1):
            for ipb_sign in (-1, 1):
                stress = combined_stress_MPa(
                    tension,
                    opb_moment,
                    ipb_moment,
                    config,
                    hotspot=hotspot,
                    opb_sign=opb_sign,
                    ipb_sign=ipb_sign,
                )
                damage, cycles = fatigue_damage_from_stress(stress, sn_curve, exposure_scale)
                sign_results.append(
                    {
                        "opb_sign": opb_sign,
                        "ipb_sign": ipb_sign,
                        "damage": damage,
                        "cycle_count": float(cycles["count"].sum()) if not cycles.empty else 0.0,
                        "stress_range_max_MPa": float(cycles["range"].max()) if not cycles.empty else 0.0,
                    }
                )

        governing = max(sign_results, key=lambda item: item["damage"])
        results["hotspots"][hotspot] = {
            "governing": governing,
            "all_sign_combinations": sign_results,
        }

    return results


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="BV NI604-style OPB/IPB top-chain fatigue from CSV time series."
    )
    parser.add_argument("csv", help="Input CSV with tension_kN and optional bending columns")
    parser.add_argument("--oc5-line", type=int, choices=[1, 2, 3], help="Use OC5 mooring line preset")
    parser.add_argument("--diameter-mm", type=float, default=None, help="Nominal chain diameter")
    parser.add_argument("--mbl-kN", type=float, default=None, help="Minimum breaking load")
    parser.add_argument("--pretension-kN", type=float, default=None, help="Pretension")
    parser.add_argument("--gamma-tt", type=float, default=None, help="Override hotspot C mean stress factor")
    parser.add_argument("--chain-type", default="studless", choices=["studless", "studlink"])
    parser.add_argument("--dt-s", type=float, default=None, help="Sample interval in seconds")
    parser.add_argument("--time-column", default=None, help="Column containing time in seconds")
    parser.add_argument("--exposure-years", type=float, default=None, help="Scale damage to exposure years")
    parser.add_argument("--occurrence-probability", type=float, default=1.0)
    parser.add_argument("--log10-k", type=float, default=12.436)
    parser.add_argument("--sn-m", type=float, default=3.0)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    df = pd.read_csv(args.csv)

    oc5_line = None
    if args.oc5_line is not None:
        from fps_motion_to_ni604_input import oc5_mooring_line

        oc5_line = oc5_mooring_line(args.oc5_line)

    diameter_mm = args.diameter_mm
    pretension_kN = args.pretension_kN
    gamma_tt = args.gamma_tt
    if oc5_line is not None:
        diameter_mm = diameter_mm or oc5_line.chain_bar_diameter_m * 1000.0
        pretension_kN = pretension_kN or oc5_line.pretension_N / 1000.0
        gamma_tt = gamma_tt if gamma_tt is not None else 0.95

    if diameter_mm is None:
        raise SystemExit("--diameter-mm is required unless --oc5-line is used")
    if pretension_kN is None:
        raise SystemExit("--pretension-kN is required unless --oc5-line is used")
    if args.mbl_kN is None and gamma_tt is None:
        raise SystemExit("--mbl-kN or --gamma-tt is required for hotspot C")

    config = ChainConfig(
        diameter_mm=diameter_mm,
        mbl_kN=args.mbl_kN,
        pretension_kN=pretension_kN,
        chain_type=args.chain_type,
        gamma_tt_override=gamma_tt,
    )
    exposure_seconds = (
        args.exposure_years * SECONDS_PER_YEAR if args.exposure_years is not None else None
    )
    results = fatigue_from_dataframe(
        df,
        config=config,
        sn_curve=SNCurve(log10_k=args.log10_k, m=args.sn_m),
        time_column=args.time_column,
        dt_s=args.dt_s,
        exposure_seconds=exposure_seconds,
        occurrence_probability=args.occurrence_probability,
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
