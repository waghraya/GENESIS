import json
import numpy as np

import orekit_jpype as orekit
orekit.initVM()
print('Virtual Machine Initialized...')
from orekit_jpype.pyhelpers import setup_orekit_data, absolutedate_to_datetime
setup_orekit_data()
print('Orekit data identified...')
from org.orekit.attitudes import Attitude, AttitudeProvider
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions, Constants
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngleType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.propagation.integration import AdditionalDerivativesProvider
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator

from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.hipparchus.geometry.euclidean.threed import Vector3D, Rotation

from objects import Scenario, Satellite, RSO

## Setup frames and time
# J2000 frame
EME2000 = FramesFactory.getEME2000()
gcrf = FramesFactory.getGCRF()
earth_frame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
utc = TimeScalesFactory.getUTC()

# TODO: move parsing somewhere else
## --------------------------------Parse input file--------------------------------
with open('Inputs/Scenario.json', 'r', encoding='utf8') as file:
    scenario_data = json.load(file)

scenario_epoch = AbsoluteDate(scenario_data['epoch'],utc)
scenario = Scenario(
    name=scenario_data['name'],
    epoch=scenario_epoch,
    duration_days=scenario_data['duration_days'],
    time_step_s=scenario_data['time_step_s']
)

with open('Inputs/Satellite.json', 'r', encoding='utf-8') as file:
    satellite_data = json.load(file)

semi_major_axis_km = satellite_data['orbit']['elements']['semi_major_axis_km']
eccentricity = satellite_data['orbit']['elements']['eccentricity']
inclination_deg = satellite_data['orbit']['elements']['inclination_deg']
raan_degrees = satellite_data['orbit']['elements']['raan_deg']
arg_of_perigee_deg = satellite_data['orbit']['elements']['arg_of_perigee_deg']
true_anomaly_deg = satellite_data['orbit']['elements']['true_anomaly_deg']
satellite_epoch = AbsoluteDate(satellite_data['epoch'],utc)
initial_orbit = KeplerianOrbit(
    semi_major_axis_km * 1000.0,    # semi-major axis in m
    eccentricity,                   # eccentricity
    np.deg2rad(inclination_deg),    # inclination
    np.deg2rad(arg_of_perigee_deg), # argument of perigee
    np.deg2rad(raan_degrees),       # raan
    np.deg2rad(true_anomaly_deg),   # true anomaly
    PositionAngleType.TRUE,         # tells which anomaly type ^
    FramesFactory.getGCRF(),        # frame
    satellite_epoch,                # epoch
    Constants.WGS84_EARTH_MU        # mu
)
sat = Satellite(
    name=satellite_data['satellite']['name'],
    id=int(satellite_data['satellite']['id']),
    international_designator=satellite_data['satellite']['international_designator'],
    epoch=satellite_epoch,
    orbit=initial_orbit,
    mass_kg=satellite_data['physical_properties']['mass_kg'],
    drag_coefficient=satellite_data['physical_properties']['drag_coefficient'],
    cross_sectional_area_m2=satellite_data['physical_properties']['cross_sectional_area_m2'],
    srp_coefficient=satellite_data['physical_properties']['srp_coefficient'],
    srp_area_m2=satellite_data['physical_properties']['srp_area_m2'],
    force_models=satellite_data['propagator']['force_models']
)

with open('Inputs/RSO.json', 'r', encoding='utf8') as file:
    rso_data = json.load(file)
rso_list = []
rso_propagators = []
for rso in rso_data['rsos']:
     rso_obj = RSO(name=rso['name'], tle=TLE(rso['tle_line1'], rso['tle_line2']))
     rso_list.append(rso_obj)
     rso_propagators.append(TLEPropagator.selectExtrapolator(rso_obj.tle))

## Integrator setup parameters
min_step = 1e-3 * scenario.time_step_s
max_step = scenario.time_step_s
pos_tolerance = 1.0

tolerances = NumericalPropagator.tolerances(pos_tolerance, sat.orbit, sat.orbit.getType())
integrator = DormandPrince853Integrator(
    min_step, max_step,
    tolerances[0],  # absolute tolerances
    tolerances[1]   # relative tolerances
)

## --------------------------------Propagator Setup--------------------------------
initial_state = SpacecraftState(sat.orbit, sat.mass_kg)
propagator = NumericalPropagator(integrator)
propagator.setOrbitType(OrbitType.CARTESIAN)
propagator.setInitialState(initial_state)

# Add an Earth gravity field force model (degree/order 8x8, for example)
earth = OneAxisEllipsoid(
    Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
    Constants.WGS84_EARTH_FLATTENING,
    earth_frame
)
gravity_provider = GravityFieldFactory.getNormalizedProvider(8, 8)
gravity_force = HolmesFeatherstoneAttractionModel(earth_frame, gravity_provider)
propagator.addForceModel(gravity_force)

## --------------------------------Propagate--------------------------------
current_date = sat.epoch
end_date = sat.epoch.shiftedBy(scenario.duration_days * 86400.0)

t = 0.0
while t < scenario.duration_days * 86400.0:
    current_date = current_date.shiftedBy(scenario.time_step_s)
    # Satellite Propagation
    satellite_state = propagator.propagate(current_date)
    satellite_pv = satellite_state.getPVCoordinates()
    # RSO Propagation
    rso_pvs = {}
    for rso, rso_prop in zip(rso_list, rso_propagators):
        rso_state = rso_prop.propagate(current_date)
        rso_pvs[rso.name] = rso_state.getPVCoordinates(gcrf)

    # Attitude Dynamics
    los = rso_pvs[rso_list[0].name].getPosition().subtract(satellite_pv.getPosition()).normalize()
    desired_rotation = Rotation(Vector3D.PLUS_K, los)
    # this is a really bad implementation that reconstructs state, need to implement controls later and propagate the input
    new_state = SpacecraftState(satellite_state.getOrbit(),
                                #eventually fill Vector3D.ZERO with attitude control outputs
                                Attitude(current_date, gcrf, desired_rotation, Vector3D.ZERO, Vector3D.ZERO)
                                )
    propagator.resetInitialState(new_state)


    #print("Date:", absolutedate_to_datetime(state.getDate()))
    #print("Position (m):", satellite_pv.getPosition())
    #print("Velocity (m/s):", satellite_pv.getVelocity())
    print(t)
    t += scenario.time_step_s
