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

def propagate(frames: dict, scenario: Scenario, sat: Satellite, rso_list: list[RSO], rso_propagators):
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
        frames['earth_frame']
    )
    gravity_provider = GravityFieldFactory.getNormalizedProvider(8, 8)
    gravity_force = HolmesFeatherstoneAttractionModel(frames['earth_frame'], gravity_provider)
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
            rso_pvs[rso.name] = rso_state.getPVCoordinates(frames['gcrf'])

        # Attitude Dynamics
        los = rso_pvs[rso_list[0].name].getPosition().subtract(satellite_pv.getPosition()).normalize()
        desired_rotation = Rotation(Vector3D.PLUS_K, los)
        # this is a really bad implementation that reconstructs state, need to implement controls later and propagate the input
        new_state = SpacecraftState(satellite_state.getOrbit(),
                                    #eventually fill Vector3D.ZERO with attitude control outputs
                                    Attitude(current_date, frames['gcrf'], desired_rotation, Vector3D.ZERO, Vector3D.ZERO)
                                    )
        propagator.resetInitialState(new_state)


        #print("Date:", absolutedate_to_datetime(state.getDate()))
        #print("Position (m):", satellite_pv.getPosition())
        #print("Velocity (m/s):", satellite_pv.getVelocity())
        print(t)
        t += scenario.time_step_s

    return 0
