"""
Convert FPSO 6DOF motion time series to screening-level NI604 OPB/IPB input.

This helper estimates fairlead relative angles from FPSO motions and maps
them to top-chain interlink angles. It is suitable for early screening and
workflow checks. For class/design use, replace the angle transfer model with
project-specific FCS tests, top-chain beam/FEM results, or mooring analysis
post-processing.

Coordinate convention:
    x: forward
    y: port
    z: up

Expected motion columns:
    time_s optional
    surge_m, sway_m, heave_m optional, default 0
    roll_deg, pitch_deg, yaw_deg optional, default 0
    tension_kN optional. If absent, tension is estimated from line stretch.

Output columns usable by ni604_opb_fatigue.py:
    tension_kN
    opb_interlink_angle_deg
    ipb_interlink_angle_deg
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from top_chain_beam_fem import TopChainBeamConfig, solve_interlink_angles_deg


@dataclass(frozen=True)
class MooringGeometry:
    """Geometry required for a single mooring-line fairlead."""

    fairlead_body_m: tuple[float, float, float]
    anchor_global_m: tuple[float, float, float]
    vessel_origin_global_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    top_link_index: int = 1
    top_link_count: int = 20
    angle_decay_ratio: float = 0.75

    def interlink_angle_factor(self) -> float:
        """Share of fairlead relative angle assigned to the selected interlink."""

        if self.top_link_index < 1:
            raise ValueError("top_link_index must be 1 or greater")
        if self.top_link_count < self.top_link_index:
            raise ValueError("top_link_count must be >= top_link_index")
        if not 0.0 < self.angle_decay_ratio < 1.0:
            raise ValueError("angle_decay_ratio must be between 0 and 1")

        r = self.angle_decay_ratio
        denominator = 1.0 - r**self.top_link_count
        return (1.0 - r) * r ** (self.top_link_index - 1) / denominator


@dataclass(frozen=True)
class TensionEstimate:
    """Optional quasi-static tension estimate when no tension column exists."""

    pretension_kN: float
    line_stiffness_kN_per_m: float
    minimum_tension_kN: float = 0.0


@dataclass(frozen=True)
class OC5MooringLine:
    """OC5-DeepCwind mooring properties for one line."""

    line_id: int
    name: str
    heading_deg: float
    fairlead_radius_m: float
    fairlead_depth_m: float
    anchor_radius_m: float
    anchor_depth_m: float
    unstretched_length_m: float
    link_inner_length_m: float
    chain_bar_diameter_m: float
    volume_equivalent_diameter_m: float
    mass_density_kg_per_m: float
    submerged_weight_N: float
    extensional_stiffness_N: float
    pretension_N: float

    @property
    def geometry(self) -> MooringGeometry:
        angle = math.radians(self.heading_deg)
        # OC5 line headings are represented as radial fairlead/anchor
        # locations around the platform center in the body/global XY plane.
        fairlead = (
            self.fairlead_radius_m * math.cos(angle),
            self.fairlead_radius_m * math.sin(angle),
            -self.fairlead_depth_m,
        )
        anchor = (
            self.anchor_radius_m * math.cos(angle),
            self.anchor_radius_m * math.sin(angle),
            -self.anchor_depth_m,
        )
        return MooringGeometry(fairlead_body_m=fairlead, anchor_global_m=anchor)

    @property
    def tension_estimate(self) -> TensionEstimate:
        return TensionEstimate(
            pretension_kN=self.pretension_N / 1000.0,
            line_stiffness_kN_per_m=(
                self.extensional_stiffness_N / self.unstretched_length_m / 1000.0
            ),
            minimum_tension_kN=0.0,
        )


OC5_MOORING_LINES = {
    1: OC5MooringLine(
        line_id=1,
        name="SBA",
        heading_deg=0.0,
        fairlead_radius_m=40.868,
        fairlead_depth_m=14.0,
        anchor_radius_m=837.6,
        anchor_depth_m=200.0,
        unstretched_length_m=835.5,
        link_inner_length_m=0.4203,
        chain_bar_diameter_m=0.0779,
        volume_equivalent_diameter_m=0.1369,
        mass_density_kg_per_m=125.6,
        submerged_weight_N=9.051e5,
        extensional_stiffness_N=7.520e8,
        pretension_N=1.107e6,
    ),
    2: OC5MooringLine(
        line_id=2,
        name="BOW",
        heading_deg=120.0,
        fairlead_radius_m=40.868,
        fairlead_depth_m=14.0,
        anchor_radius_m=837.6,
        anchor_depth_m=200.0,
        unstretched_length_m=835.5,
        link_inner_length_m=0.4203,
        chain_bar_diameter_m=0.0779,
        volume_equivalent_diameter_m=0.1398,
        mass_density_kg_per_m=125.8,
        submerged_weight_N=9.017e5,
        extensional_stiffness_N=7.461e8,
        pretension_N=1.112e6,
    ),
    3: OC5MooringLine(
        line_id=3,
        name="PSA",
        heading_deg=240.0,
        fairlead_radius_m=40.868,
        fairlead_depth_m=14.0,
        anchor_radius_m=837.6,
        anchor_depth_m=200.0,
        unstretched_length_m=835.5,
        link_inner_length_m=0.4203,
        chain_bar_diameter_m=0.0779,
        volume_equivalent_diameter_m=0.1393,
        mass_density_kg_per_m=125.4,
        submerged_weight_N=8.997e5,
        extensional_stiffness_N=7.478e8,
        pretension_N=1.148e6,
    ),
}


def oc5_mooring_line(line_id: int) -> OC5MooringLine:
    """Return full-scale OC5-DeepCwind mooring preset for line 1, 2, or 3."""

    if line_id not in OC5_MOORING_LINES:
        raise ValueError("OC5 line_id must be 1, 2, or 3")
    return OC5_MOORING_LINES[line_id]


def _column_or_zero(df: pd.DataFrame, column: str) -> np.ndarray:
    if column in df.columns:
        return df[column].to_numpy(dtype=float)
    return np.zeros(len(df), dtype=float)


def _unit(vector: Iterable[float] | np.ndarray) -> np.ndarray:
    arr = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(arr)
    if norm <= 0.0:
        raise ValueError("Cannot normalize a zero vector")
    return arr / norm


def _rotation_matrix_body_to_global(roll_deg: float, pitch_deg: float, yaw_deg: float) -> np.ndarray:
    """Body-to-global rotation using Rz(yaw) * Ry(pitch) * Rx(roll)."""

    roll = math.radians(roll_deg)
    pitch = math.radians(pitch_deg)
    yaw = math.radians(yaw_deg)

    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)

    rx = np.array([[1.0, 0.0, 0.0], [0.0, cr, -sr], [0.0, sr, cr]])
    ry = np.array([[cp, 0.0, sp], [0.0, 1.0, 0.0], [-sp, 0.0, cp]])
    rz = np.array([[cy, -sy, 0.0], [sy, cy, 0.0], [0.0, 0.0, 1.0]])
    return rz @ ry @ rx


def _reference_chain_basis(
    fairlead_body_m: np.ndarray,
    anchor_global_m: np.ndarray,
    origin_global_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return nominal line, in-plane, and out-of-plane unit vectors in body axes."""

    line_axis = _unit(anchor_global_m - (origin_global_m + fairlead_body_m))
    global_up = np.array([0.0, 0.0, 1.0])

    # The vertical plane containing the nominal line defines "in-plane".
    # Its normal direction defines "out-of-plane" for the fairlead angle.
    out_of_plane = np.cross(global_up, line_axis)
    if np.linalg.norm(out_of_plane) < 1.0e-9:
        out_of_plane = np.array([0.0, 1.0, 0.0])
    out_of_plane = _unit(out_of_plane)
    in_plane = _unit(np.cross(out_of_plane, line_axis))

    return line_axis, in_plane, out_of_plane


def fps_motion_to_fairlead_angles(
    df: pd.DataFrame,
    geometry: MooringGeometry,
) -> pd.DataFrame:
    """Compute fairlead relative in-plane and out-of-plane angles."""

    fairlead_body = np.asarray(geometry.fairlead_body_m, dtype=float)
    anchor_global = np.asarray(geometry.anchor_global_m, dtype=float)
    origin_global = np.asarray(geometry.vessel_origin_global_m, dtype=float)
    line_axis_body, in_plane_body, out_of_plane_body = _reference_chain_basis(
        fairlead_body,
        anchor_global,
        origin_global,
    )

    surge = _column_or_zero(df, "surge_m")
    sway = _column_or_zero(df, "sway_m")
    heave = _column_or_zero(df, "heave_m")
    roll = _column_or_zero(df, "roll_deg")
    pitch = _column_or_zero(df, "pitch_deg")
    yaw = _column_or_zero(df, "yaw_deg")

    rows: list[dict[str, float]] = []
    for i in range(len(df)):
        rotation = _rotation_matrix_body_to_global(roll[i], pitch[i], yaw[i])
        vessel_translation = origin_global + np.array([surge[i], sway[i], heave[i]])
        fairlead_global = vessel_translation + rotation @ fairlead_body

        # The actual top-line direction points from the moved fairlead toward
        # the fixed anchor. The nominal chain-stopper basis rotates with FPSO.
        actual_line_axis = _unit(anchor_global - fairlead_global)
        fcs_line_axis = rotation @ line_axis_body
        fcs_in_plane = rotation @ in_plane_body
        fcs_out_of_plane = rotation @ out_of_plane_body

        axial_component = float(np.dot(actual_line_axis, fcs_line_axis))
        in_plane_component = float(np.dot(actual_line_axis, fcs_in_plane))
        out_of_plane_component = float(np.dot(actual_line_axis, fcs_out_of_plane))

        rows.append(
            {
                "fairlead_x_m": float(fairlead_global[0]),
                "fairlead_y_m": float(fairlead_global[1]),
                "fairlead_z_m": float(fairlead_global[2]),
                "line_length_m": float(np.linalg.norm(anchor_global - fairlead_global)),
                "fairlead_angle_in_plane_deg": math.degrees(
                    # atan2 keeps the sign and stays stable for small angles.
                    math.atan2(in_plane_component, axial_component)
                ),
                "fairlead_angle_out_of_plane_deg": math.degrees(
                    math.atan2(out_of_plane_component, axial_component)
                ),
            }
        )

    return pd.DataFrame(rows)


def build_ni604_input_from_fps_motion(
    df: pd.DataFrame,
    geometry: MooringGeometry,
    tension_estimate: TensionEstimate | None = None,
    angle_model: str = "decay",
    beam_config: TopChainBeamConfig | None = None,
    include_all_interlinks: bool = False,
) -> pd.DataFrame:
    """Create NI604 fatigue input columns from FPSO motion time series."""

    angles = fps_motion_to_fairlead_angles(df, geometry)
    output = pd.DataFrame(index=df.index)

    if "time_s" in df.columns:
        output["time_s"] = df["time_s"].to_numpy(dtype=float)

    if "tension_kN" in df.columns:
        output["tension_kN"] = df["tension_kN"].to_numpy(dtype=float)
    elif tension_estimate is not None:
        reference_length = float(angles["line_length_m"].iloc[0])
        extension = angles["line_length_m"].to_numpy(dtype=float) - reference_length
        # This is a quasi-static screening estimate. A production workflow
        # should use mooring-analysis tension if it is available.
        tension = (
            tension_estimate.pretension_kN
            + tension_estimate.line_stiffness_kN_per_m * extension
        )
        output["tension_kN"] = np.maximum(tension, tension_estimate.minimum_tension_kN)
    else:
        raise ValueError(
            "Input needs tension_kN or a TensionEstimate for screening tension calculation"
        )

    if angle_model == "decay":
        # Fast screening model: prescribe how much of the fairlead angle is
        # assigned to the selected top-chain interlink.
        factor = geometry.interlink_angle_factor()
        output["opb_interlink_angle_deg"] = angles["fairlead_angle_out_of_plane_deg"] * factor
        output["ipb_interlink_angle_deg"] = angles["fairlead_angle_in_plane_deg"] * factor
    elif angle_model == "fem":
        if beam_config is None:
            beam_config = TopChainBeamConfig(link_count=geometry.top_link_count)
        if geometry.top_link_index >= beam_config.link_count:
            raise ValueError(
                "top_link_index selects an interlink and must be less than FEM link_count"
            )
        # More physical screening model: solve a 1D beam/rotational-spring
        # top-chain model in the OPB and IPB planes independently.
        opb_interlinks = solve_interlink_angles_deg(
            angles["fairlead_angle_out_of_plane_deg"].to_numpy(dtype=float),
            output["tension_kN"].to_numpy(dtype=float),
            beam_config,
        )
        ipb_interlinks = solve_interlink_angles_deg(
            angles["fairlead_angle_in_plane_deg"].to_numpy(dtype=float),
            output["tension_kN"].to_numpy(dtype=float),
            beam_config,
        )
        selected = geometry.top_link_index - 1
        output["opb_interlink_angle_deg"] = opb_interlinks[:, selected]
        output["ipb_interlink_angle_deg"] = ipb_interlinks[:, selected]
        if include_all_interlinks:
            for interlink in range(beam_config.link_count - 1):
                output[f"opb_interlink_{interlink + 1:02d}_deg"] = opb_interlinks[:, interlink]
                output[f"ipb_interlink_{interlink + 1:02d}_deg"] = ipb_interlinks[:, interlink]
    else:
        raise ValueError("angle_model must be 'decay' or 'fem'")

    return pd.concat([output, angles], axis=1)


def make_synthetic_fps_motion(duration_s: float = 3600.0, dt_s: float = 0.5) -> pd.DataFrame:
    """Generate a small synthetic FPSO motion/tension time series for examples."""

    oc5_line_1 = oc5_mooring_line(1)
    pretension_kN = oc5_line_1.pretension_N / 1000.0
    time = np.arange(0.0, duration_s + 0.5 * dt_s, dt_s)
    wave = 2.0 * math.pi * time / 12.5
    slow = 2.0 * math.pi * time / 120.0

    return pd.DataFrame(
        {
            "time_s": time,
            "surge_m": 4.0 * np.sin(slow) + 0.7 * np.sin(wave + 0.2),
            "sway_m": 1.2 * np.sin(0.8 * slow + 0.7) + 0.4 * np.sin(0.9 * wave),
            "heave_m": 0.8 * np.sin(wave - 0.4),
            "roll_deg": 1.2 * np.sin(0.95 * wave + 0.3),
            "pitch_deg": 1.8 * np.sin(1.05 * wave - 0.5),
            "yaw_deg": 2.5 * np.sin(0.7 * slow + 0.1),
            "tension_kN": (
                pretension_kN + 90.0 * np.sin(slow + 0.4) + 35.0 * np.sin(wave)
            ),
        }
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create NI604 OPB/IPB angle input from FPSO 6DOF motion CSV."
    )
    parser.add_argument("input_csv", nargs="?", help="FPSO motion CSV")
    parser.add_argument("output_csv", nargs="?", help="Output CSV for ni604_opb_fatigue.py")
    parser.add_argument("--make-demo", action="store_true", help="Write synthetic demo CSV files")
    parser.add_argument("--demo-duration-s", type=float, default=3600.0)
    parser.add_argument("--demo-dt-s", type=float, default=0.5)
    parser.add_argument("--oc5-line", type=int, choices=[1, 2, 3], help="Use OC5 mooring line preset")
    parser.add_argument("--fairlead", nargs=3, type=float, metavar=("X", "Y", "Z"))
    parser.add_argument("--anchor", nargs=3, type=float, metavar=("X", "Y", "Z"))
    parser.add_argument("--origin", nargs=3, type=float, default=(0.0, 0.0, 0.0))
    parser.add_argument("--top-link-index", type=int, default=1)
    parser.add_argument("--top-link-count", type=int, default=20)
    parser.add_argument("--angle-decay-ratio", type=float, default=0.75)
    parser.add_argument("--angle-model", default="decay", choices=["decay", "fem"])
    parser.add_argument("--write-all-interlinks", action="store_true")
    parser.add_argument("--beam-link-length-m", type=float, default=None)
    parser.add_argument("--beam-bar-diameter-m", type=float, default=None)
    parser.add_argument("--beam-e-pa", type=float, default=None)
    parser.add_argument("--beam-end-rotation", default="pinned", choices=["pinned", "aligned"])
    parser.add_argument("--beam-max-iterations", type=int, default=8)
    parser.add_argument("--beam-convergence-tol-deg", type=float, default=1.0e-5)
    parser.add_argument("--pretension-kN", type=float, default=None)
    parser.add_argument("--line-stiffness-kN-per-m", type=float, default=0.0)
    parser.add_argument("--minimum-tension-kN", type=float, default=0.0)
    return parser


def _beam_config_from_args(
    args: argparse.Namespace,
    oc5_line: OC5MooringLine | None,
) -> TopChainBeamConfig:
    link_length_m = args.beam_link_length_m
    bar_diameter_m = args.beam_bar_diameter_m
    elastic_modulus_pa = args.beam_e_pa

    if oc5_line is not None:
        # OC5 uses a brass chain model; keep E overridable because real steel
        # mooring chain should use a different elastic modulus.
        link_length_m = link_length_m or oc5_line.link_inner_length_m
        bar_diameter_m = bar_diameter_m or oc5_line.chain_bar_diameter_m
        elastic_modulus_pa = elastic_modulus_pa or 100.0e9
    else:
        if bar_diameter_m is None:
            bar_diameter_m = 0.12
        link_length_m = link_length_m or 6.0 * bar_diameter_m
        elastic_modulus_pa = elastic_modulus_pa or 210.0e9

    return TopChainBeamConfig(
        link_count=args.top_link_count,
        link_length_m=link_length_m,
        bar_diameter_m=bar_diameter_m,
        elastic_modulus_pa=elastic_modulus_pa,
        end_rotation=args.beam_end_rotation,
        max_iterations=args.beam_max_iterations,
        convergence_tol_deg=args.beam_convergence_tol_deg,
    )


def main() -> None:
    args = _build_parser().parse_args()

    if args.make_demo:
        demo_motion = make_synthetic_fps_motion(
            duration_s=args.demo_duration_s,
            dt_s=args.demo_dt_s,
        )
        motion_path = Path("fps_motion_demo.csv")
        output_path = Path(
            "ni604_input_demo.csv"
            if args.angle_model == "decay"
            else f"ni604_input_demo_{args.angle_model}.csv"
        )
        oc5_line = oc5_mooring_line(args.oc5_line or 1)
        geometry = MooringGeometry(
            fairlead_body_m=oc5_line.geometry.fairlead_body_m,
            anchor_global_m=oc5_line.geometry.anchor_global_m,
            top_link_index=args.top_link_index,
            top_link_count=args.top_link_count,
            angle_decay_ratio=args.angle_decay_ratio,
        )
        beam_config = (
            _beam_config_from_args(args, oc5_line) if args.angle_model == "fem" else None
        )
        ni604_input = build_ni604_input_from_fps_motion(
            demo_motion,
            geometry,
            angle_model=args.angle_model,
            beam_config=beam_config,
            include_all_interlinks=args.write_all_interlinks,
        )
        demo_motion.to_csv(motion_path, index=False)
        ni604_input.to_csv(output_path, index=False)
        print(f"Wrote {motion_path} and {output_path}")
        return

    if not args.input_csv or not args.output_csv:
        raise SystemExit("input_csv and output_csv are required unless --make-demo is used")
    oc5_line = oc5_mooring_line(args.oc5_line) if args.oc5_line else None
    if oc5_line is None and (args.fairlead is None or args.anchor is None):
        raise SystemExit("--fairlead and --anchor are required")

    df = pd.read_csv(args.input_csv)
    base_geometry = oc5_line.geometry if oc5_line is not None else MooringGeometry(
        fairlead_body_m=tuple(args.fairlead),
        anchor_global_m=tuple(args.anchor),
    )
    geometry = MooringGeometry(
        fairlead_body_m=base_geometry.fairlead_body_m,
        anchor_global_m=base_geometry.anchor_global_m,
        vessel_origin_global_m=tuple(args.origin),
        top_link_index=args.top_link_index,
        top_link_count=args.top_link_count,
        angle_decay_ratio=args.angle_decay_ratio,
    )
    tension_estimate = None
    if "tension_kN" not in df.columns:
        if oc5_line is not None and args.pretension_kN is None:
            tension_estimate = oc5_line.tension_estimate
        elif args.pretension_kN is None:
            raise SystemExit("--pretension-kN is required when input has no tension_kN column")
        else:
            tension_estimate = TensionEstimate(
                pretension_kN=args.pretension_kN,
                line_stiffness_kN_per_m=args.line_stiffness_kN_per_m,
                minimum_tension_kN=args.minimum_tension_kN,
            )

    beam_config = _beam_config_from_args(args, oc5_line) if args.angle_model == "fem" else None
    output = build_ni604_input_from_fps_motion(
        df,
        geometry,
        tension_estimate,
        angle_model=args.angle_model,
        beam_config=beam_config,
        include_all_interlinks=args.write_all_interlinks,
    )
    output.to_csv(args.output_csv, index=False)
    print(f"Wrote {args.output_csv}")


if __name__ == "__main__":
    main()
