
"""
cepa_surface_loading.py
Minimal calculateur CEPA "Surface Loading" – version MVP.
ATTENTION: Certains paramètres sol (Kb, Kz) sont laissés en entrée explicite jusqu’à calibration.
"""

from dataclasses import dataclass
from typing import Literal, Optional
import math

Units = Literal["SI", "IMP"]

@dataclass
class Pipe:
    D: float          # outside diameter [m (SI) / in (IMP)]
    t: float          # wall thickness [m / in]
    MOP: float        # max operating pressure [Pa / psig]
    SMYS: float       # yield strength [Pa / psi]
    dT: float         # temperature differential [°C / °F]
    units: Units = "SI"

@dataclass
class Soil:
    rho: float        # density [kg/m^3 / lb/ft^3]
    H: float          # cover depth [m / ft]
    theta_deg: float  # bedding angle [deg]
    Eprime: float     # modulus of soil reaction E' [Pa / psi]
    phi_deg: Optional[float] = None  # friction angle, needed for Trap Door if used
    load_equation: Literal["prism", "trap_door"] = "prism"
    units: Units = "SI"

@dataclass
class VehicleWheels:
    # Two- or three-axle wheel vehicle; provide arrays matching number of axles
    axle_loads: list[float]      # per-axle gross load [N / lb]
    axle_separations: list[float]# distances between consecutive axles [m / ft]; len = n_axles-1
    axle_width: float            # lateral spacing between left/right tire centerlines [m / ft]
    contact_widths: list[float]  # tire contact widths per axle [m / in]; duals treated as one
    tire_pressures: list[float]  # per-axle tire pressure [Pa / psi]
    units: Units = "SI"

@dataclass
class VehicleTracks:
    gross_load: float    # total vehicle operating weight [N / lb]
    track_length: float  # ground contact length [m / ft]
    track_sep: float     # distance between track centerlines [m / ft]
    track_width: float   # width of one track [m / in]
    units: Units = "SI"

@dataclass
class Misc:
    vehicle_type: Literal["highway", "farm_construction", "track"] = "highway"
    pavement: Literal["rigid", "flexible"] = "flexible"
    eqv: Literal["tresca", "von_mises"] = "tresca"
    poisson: float = 0.3
    alpha: float = 6.67e-6  # thermal expansion [1/°C]
    impact_factor: Optional[float] = None  # if None, pick by vehicle_type+pavement

# --- Unit helpers ---
def to_si(value: float, src: Units, what: str) -> float:
    """Convert basic inputs to SI for internal math."""
    if src == "SI":
        return value
    # IMP -> SI conversions
    if what in ("length_m",):
        return value * 0.3048  # ft->m
    if what in ("length_in",):
        return value * 0.0254   # in->m
    if what in ("pressure_psig",):
        # psig -> Pa (gauge); we'll treat as Pa for stress deltas
        return value * 6894.75729
    if what in ("pressure_psi",):
        return value * 6894.75729
    if what in ("stress_psi",):
        return value * 6894.75729
    if what in ("density_lbft3",):
        return value * 16.018463  # lb/ft3 -> kg/m3
    if what in ("force_lb",):
        return value * 4.4482216152605  # lb(force) -> N
    if what in ("tempF",):
        return (value - 32) * 5.0/9.0
    raise ValueError(f"Unknown conversion key {what}")

def impact_factor(m: Misc) -> float:
    if m.impact_factor is not None:
        return m.impact_factor
    # Heuristics inspired by manual notes (order of magnitude; refine per project policy)
    if m.vehicle_type == "highway":
        return 1.15 if m.pavement == "flexible" else 1.00
    if m.vehicle_type == "farm_construction":
        return 1.20 if m.pavement == "flexible" else 1.05
    if m.vehicle_type == "track":
        return 1.25 if m.pavement == "flexible" else 1.10
    return 1.15

# --- Core equations (SI internal) ---

def boussinesq_pressure(F_point: float, H: float, D: float, d_offset: float = 0.0) -> float:
    """Pressure at pipe crown from a point load at surface (Boussinesq). Eq (1).
    liveP = (W_live / D), with W_live from point load F. We implement the direct form.

    Units: F [N], H [m], D [m] -> P [Pa].
    """
    denom = (H**2 + d_offset**2) ** 2.5
    return (F_point / (2.0 * math.pi)) * (3.0 * H**3 / denom) / D

def soil_prism_pressure(rho: float, H: float, D: float) -> float:
    """Soil prism load pressure on pipe crown. Eq (2).
    In SI: P_soil = rho * g * H.
    """
    g = 9.80665
    return rho * g * H

def hoop_from_internal_pressure(P_internal: float, D: float, t: float) -> float:
    """Hoop stress due to internal pressure. Eq (5): sigma_H = P*D/(2t)."""
    return P_internal * D / (2.0 * t)

def longitudinal_from_pressure(hoop_pressure_stress: float, nu: float = 0.3) -> float:
    """Longitudinal stress due to internal pressure (Poisson effect). Eq (6)."""
    return nu * hoop_pressure_stress

def longitudinal_from_soil(hoop_soil_stress: float, nu: float = 0.3) -> float:
    """Longitudinal stress due to soil load (Poisson effect). Eq (7)."""
    return nu * hoop_soil_stress

def cepa_hoop_bending(P_live: float, D: float, t: float, E: float, Eprime: float, P_internal: float, Kb: float, Kz: float) -> float:
    """CEPA hoop bending component (Eq (3) / (4)).
    sigma_H_bending = (P * (D/t)^3) * [ Kb*(E/E') + 0.915*(P_internal*D/(E*t)) + Kz*(E/E')*(P_internal*D/(E*t)) ]
    NOTE: Kb, Kz must be calibrated from project guidance.
    """
    Dt = D / t
    term = (Kb * (E / Eprime)) + 0.915 * (P_internal * D / (E * t)) + (Kz * (E / Eprime) * (P_internal * D / (E * t)))
    return P_live * (Dt ** 3) * term

def hetenyi_lambda(E: float, I: float, D: float, Eprime: float, theta_deg: float) -> float:
    """Characteristic length lambda for beam on elastic foundation (Eq (11))."""
    theta = math.radians(theta_deg)
    lam4 = (D * Eprime * theta) / (E * I)
    return lam4 ** 0.25

def longitudinal_bending_from_beam(P_pipe: float, D: float, t: float, E: float, Eprime: float, theta_deg: float) -> float:
    """Max longitudinal bending stress from beam on elastic foundation (Eqs (10)-(12))."""
    I = math.pi/64.0 * ( (D**4) - ( (D-2*t)**4 ) )
    lam = hetenyi_lambda(E, I, D, Eprime, theta_deg)
    M = (D * P_pipe) / (2.0 * (lam**2))
    c = D/2.0
    sigma = M * c / I
    return sigma

def thermal_longitudinal(alpha: float, E: float, dT: float) -> float:
    return alpha * E * dT

def equivalent_stress(sigma_H: float, sigma_L: float, method: Literal["tresca","von_mises"]="tresca") -> float:
    if method == "tresca":
        return max(abs(sigma_H - sigma_L), abs(sigma_H), abs(sigma_L))
    return math.sqrt(sigma_H**2 - sigma_H*sigma_L + sigma_L**2)

def example_compute_over_crown(pipe: Pipe, soil: Soil, misc: Misc, F_point_surface: float, Kb: float=0.031, Kz: float=0.0915):
    """Compute key stresses for a single surface point load directly above the crown (d=0).
    WARNING: Simplified demonstrator. Real footprints require superposition of multiple points (wheels/tracks/grid).
    """
    # Convert inputs to SI
    if pipe.units == "IMP":
        D = pipe.D * 0.0254
        t = pipe.t * 0.0254
        Pmop = pipe.MOP * 6894.75729
        SMYS = pipe.SMYS * 6894.75729
        dT = (pipe.dT - 32) * 5.0/9.0
    else:
        D, t, Pmop, SMYS, dT = pipe.D, pipe.t, pipe.MOP, pipe.SMYS, pipe.dT

    if soil.units == "IMP":
        rho = soil.rho * 16.018463
        H = soil.H * 0.3048
        Eprime = soil.Eprime * 6894.75729
    else:
        rho, H, Eprime = soil.rho, soil.H, soil.Eprime

    if pipe.units == "IMP":
        F = F_point_surface * 4.4482216152605
    else:
        F = F_point_surface

    # Constants
    E = 200e9  # Pa (steel)
    nu = misc.poisson
    IF = impact_factor(misc)

    # Pressures
    P_live = boussinesq_pressure(F_point=F*IF, H=H, D=D, d_offset=0.0)
    P_soil = soil_prism_pressure(rho=rho, H=H, D=D)

    # Hoop components
    sigma_H_int = hoop_from_internal_pressure(Pmop, D, t)
    sigma_H_live_bend = cepa_hoop_bending(P_live, D, t, E, Eprime, Pmop, Kb, Kz)
    sigma_H_soil_bend = cepa_hoop_bending(P_soil, D, t, E, Eprime, 0.0, Kb, Kz)

    sigma_H_total = sigma_H_int + sigma_H_live_bend + sigma_H_soil_bend

    # Longitudinal components
    sigma_L_from_int = longitudinal_from_pressure(sigma_H_int, nu)
    sigma_L_from_soil = longitudinal_from_soil(sigma_H_soil_bend, nu)
    sigma_L_beam = longitudinal_bending_from_beam(P_pipe=P_live, D=D, t=t, E=E, Eprime=Eprime, theta_deg=soil.theta_deg)
    sigma_L_thermal = thermal_longitudinal(misc.alpha, E, dT)
    sigma_L_total = sigma_L_from_int + sigma_L_from_soil + sigma_L_beam + sigma_L_thermal

    sigma_eq = equivalent_stress(sigma_H_total, sigma_L_total, misc.eqv)

    return {
        "P_live_Pa": P_live,
        "P_soil_Pa": P_soil,
        "sigma_H_int_Pa": sigma_H_int,
        "sigma_H_live_bend_Pa": sigma_H_live_bend,
        "sigma_H_soil_bend_Pa": sigma_H_soil_bend,
        "sigma_H_total_Pa": sigma_H_total,
        "sigma_L_beam_Pa": sigma_L_beam,
        "sigma_L_total_Pa": sigma_L_total,
        "sigma_eq_Pa": sigma_eq,
        "SMYS_Pa": SMYS,
        "utilization_eq": sigma_eq / SMYS if SMYS>0 else None,
    }
