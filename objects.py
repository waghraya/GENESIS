from org.orekit.propagation.analytical.tle import TLE
from org.orekit.orbits import KeplerianOrbit
from org.orekit.time import AbsoluteDate
from dataclasses import dataclass

@dataclass
class Scenario:
    name: str
    epoch: AbsoluteDate
    duration_days: float
    time_step_s: float

@dataclass
class Satellite:
    name: str
    id: int
    international_designator: str
    epoch: AbsoluteDate
    orbit: KeplerianOrbit
    mass_kg: float
    drag_coefficient: float
    cross_sectional_area_m2: float
    srp_coefficient: float
    srp_area_m2: float
    force_models: tuple[str,...]

@dataclass(frozen=True)
class RSO:
    name: str
    tle: TLE

