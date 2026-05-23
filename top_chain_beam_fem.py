"""
1D beam / rotational-spring model for top-chain interlink angles.

The model is a screening-level local FEM:
    - each chain link is an Euler-Bernoulli beam element,
    - adjacent links are connected by nonlinear rotational springs,
    - spring stiffness is iterated from the NI604 interlink moment law,
    - OPB and IPB planes are solved independently.

It is not a substitute for project-specific chain-stopper tests or a 3D
contact FEM model, but it gives a more mechanical angle-transfer model than
a prescribed exponential decay.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np

from ni604_opb_fatigue import ni604_interlink_moment_Nmm


@dataclass(frozen=True)
class TopChainBeamConfig:
    # The default values are tuned for the OC5 brass-chain screening example.
    # Override bar diameter, link length, and E for project steel chain cases.
    link_count: int = 20
    link_length_m: float = 0.42
    bar_diameter_m: float = 0.0779
    elastic_modulus_pa: float = 100.0e9
    friction_coefficient: float = 0.30
    end_rotation: str = "pinned"
    max_iterations: int = 8
    convergence_tol_deg: float = 1.0e-5
    stiffness_relaxation: float = 0.45
    minimum_joint_stiffness_Nm_per_rad: float = 1.0e3

    @property
    def diameter_mm(self) -> float:
        return self.bar_diameter_m * 1000.0

    @property
    def flexural_rigidity_Nm2(self) -> float:
        second_moment = math.pi * self.bar_diameter_m**4 / 64.0
        return self.elastic_modulus_pa * second_moment


def _beam_bending_stiffness(ei: float, length: float) -> np.ndarray:
    l = length
    return ei / l**3 * np.array(
        [
            [12.0, 6.0 * l, -12.0, 6.0 * l],
            [6.0 * l, 4.0 * l**2, -6.0 * l, 2.0 * l**2],
            [-12.0, -6.0 * l, 12.0, -6.0 * l],
            [6.0 * l, 2.0 * l**2, -6.0 * l, 4.0 * l**2],
        ],
        dtype=float,
    )


def _beam_tension_geometric_stiffness(tension_N: float, length: float) -> np.ndarray:
    """Geometric stiffness for a beam-column under positive axial tension."""

    l = length
    return tension_N / (30.0 * l) * np.array(
        [
            [36.0, 3.0 * l, -36.0, 3.0 * l],
            [3.0 * l, 4.0 * l**2, -3.0 * l, -l**2],
            [-36.0, -3.0 * l, 36.0, -3.0 * l],
            [3.0 * l, -l**2, -3.0 * l, 4.0 * l**2],
        ],
        dtype=float,
    )


def _rotation_left_dof(link: int, node_count: int) -> int:
    return node_count + 2 * link


def _rotation_right_dof(link: int, node_count: int) -> int:
    return node_count + 2 * link + 1


def _secant_joint_stiffness_Nm_per_rad(
    angle_rad: np.ndarray,
    tension_kN: float,
    config: TopChainBeamConfig,
) -> np.ndarray:
    abs_angle_rad = np.maximum(np.abs(angle_rad), math.radians(1.0e-4))
    angle_deg = np.degrees(abs_angle_rad)
    # The nonlinear NI604 moment relation is converted to an equivalent
    # secant rotational stiffness, k = M / theta, for each interlink.
    moment_Nm = (
        ni604_interlink_moment_Nmm(
            angle_deg,
            np.full_like(angle_deg, tension_kN, dtype=float),
            diameter_mm=config.diameter_mm,
            friction_coefficient=config.friction_coefficient,
        )
        / 1000.0
    )
    stiffness = moment_Nm / abs_angle_rad
    return np.maximum(stiffness, config.minimum_joint_stiffness_Nm_per_rad)


def _assemble_stiffness(
    joint_stiffness_Nm_per_rad: np.ndarray,
    tension_kN: float,
    config: TopChainBeamConfig,
) -> np.ndarray:
    link_count = config.link_count
    node_count = link_count + 1
    dof_count = node_count + 2 * link_count

    # DOF layout:
    #   0..link_count                  : transverse node translations
    #   node_count + 2*link            : left-end rotation of each link
    #   node_count + 2*link + 1        : right-end rotation of each link
    # Rotations are duplicated at interlinks so rotational springs can carry
    # a relative angle between adjacent link ends.
    stiffness = np.zeros((dof_count, dof_count), dtype=float)
    element_k = _beam_bending_stiffness(
        config.flexural_rigidity_Nm2,
        config.link_length_m,
    )
    element_k += _beam_tension_geometric_stiffness(
        tension_kN * 1000.0,
        config.link_length_m,
    )

    for link in range(link_count):
        # Each physical link contributes beam bending plus geometric stiffness
        # from axial tension. Positive tension stiffens lateral bending.
        dofs = [
            link,
            _rotation_left_dof(link, node_count),
            link + 1,
            _rotation_right_dof(link, node_count),
        ]
        for row_local, row_global in enumerate(dofs):
            for col_local, col_global in enumerate(dofs):
                stiffness[row_global, col_global] += element_k[row_local, col_local]

    for joint in range(link_count - 1):
        # Interlink contact is represented by a rotational spring between the
        # right rotation of the previous link and the left rotation of the next.
        k = joint_stiffness_Nm_per_rad[joint]
        left = _rotation_right_dof(joint, node_count)
        right = _rotation_left_dof(joint + 1, node_count)
        stiffness[left, left] += k
        stiffness[right, right] += k
        stiffness[left, right] -= k
        stiffness[right, left] -= k

    return stiffness


def solve_interlink_angles_for_step_deg(
    fairlead_angle_deg: float,
    tension_kN: float,
    config: TopChainBeamConfig = TopChainBeamConfig(),
) -> tuple[np.ndarray, np.ndarray]:
    """Solve interlink angles and link-end rotations for one time step.

    Returns:
        interlink_angles_deg: shape (link_count - 1,)
        link_end_rotations_deg: shape (link_count, 2), left/right rotation per link
    """

    if config.link_count < 2:
        raise ValueError("link_count must be at least 2")
    if config.end_rotation not in {"pinned", "aligned"}:
        raise ValueError("end_rotation must be 'pinned' or 'aligned'")

    link_count = config.link_count
    node_count = link_count + 1
    dof_count = node_count + 2 * link_count
    fairlead_angle_rad = math.radians(fairlead_angle_deg)

    joint_angles = np.full(link_count - 1, fairlead_angle_rad / link_count, dtype=float)
    joint_stiffness = _secant_joint_stiffness_Nm_per_rad(joint_angles, tension_kN, config)
    displacements = np.zeros(dof_count, dtype=float)

    # Apply fairlead relative angle at the first link. The far node is pinned
    # in translation; optionally its last-link rotation can also be aligned.
    fixed_values = {
        0: 0.0,
        link_count: 0.0,
        _rotation_left_dof(0, node_count): fairlead_angle_rad,
    }
    if config.end_rotation == "aligned":
        fixed_values[_rotation_right_dof(link_count - 1, node_count)] = 0.0

    all_dofs = np.arange(dof_count)
    fixed = np.array(sorted(fixed_values), dtype=int)
    free = np.setdiff1d(all_dofs, fixed)
    fixed_displacements = np.array([fixed_values[dof] for dof in fixed], dtype=float)

    for _ in range(config.max_iterations):
        stiffness = _assemble_stiffness(joint_stiffness, tension_kN, config)
        rhs = -stiffness[np.ix_(free, fixed)] @ fixed_displacements
        k_free = stiffness[np.ix_(free, free)]
        try:
            displacements[free] = np.linalg.solve(k_free, rhs)
        except np.linalg.LinAlgError:
            displacements[free] = np.linalg.lstsq(k_free, rhs, rcond=None)[0]
        displacements[fixed] = fixed_displacements

        new_joint_angles = np.array(
            [
                displacements[_rotation_right_dof(joint, node_count)]
                - displacements[_rotation_left_dof(joint + 1, node_count)]
                for joint in range(link_count - 1)
            ],
            dtype=float,
        )
        max_change = float(np.max(np.abs(new_joint_angles - joint_angles)))
        joint_angles = new_joint_angles

        # Recompute nonlinear interlink stiffness from the solved angle and
        # relax the update to avoid stiffness chatter near the sliding limit.
        updated_stiffness = _secant_joint_stiffness_Nm_per_rad(joint_angles, tension_kN, config)
        joint_stiffness = (
            config.stiffness_relaxation * updated_stiffness
            + (1.0 - config.stiffness_relaxation) * joint_stiffness
        )
        if math.degrees(max_change) <= config.convergence_tol_deg:
            break

    rotations = np.array(
        [
            [
                displacements[_rotation_left_dof(link, node_count)],
                displacements[_rotation_right_dof(link, node_count)],
            ]
            for link in range(link_count)
        ],
        dtype=float,
    )

    return np.degrees(joint_angles), np.degrees(rotations)


def solve_interlink_angles_deg(
    fairlead_angle_deg: Iterable[float] | np.ndarray,
    tension_kN: Iterable[float] | np.ndarray,
    config: TopChainBeamConfig = TopChainBeamConfig(),
) -> np.ndarray:
    """Solve interlink angle time series.

    Returns an array with shape (sample_count, link_count - 1).
    """

    angles = np.asarray(fairlead_angle_deg, dtype=float)
    tensions = np.asarray(tension_kN, dtype=float)
    if angles.size != tensions.size:
        raise ValueError("fairlead_angle_deg and tension_kN must have the same length")

    output = np.zeros((angles.size, config.link_count - 1), dtype=float)
    for index, (angle, tension) in enumerate(zip(angles, tensions)):
        output[index], _ = solve_interlink_angles_for_step_deg(
            float(angle),
            float(tension),
            config,
        )
    return output
