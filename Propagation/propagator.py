## =================================================================================================
## GENESIS PROPAGATOR
## Propagates our spacecraft (orbit + attitude) and the RSOs around it, and records pointing,
## visibility and lighting at every time step.
##
## Run it with:  python run.py
## It reads the files in Inputs/ and writes the results to Outputs/<scenario name>/.
##
## Sections, in the order a run uses them:
##   1. Setup                     Orekit, file locations, reference frames, Earth/Sun/Moon
##   2. Inputs                    reads and checks Scenario.json, Satellite.json and the RSO file
##   3. Orbit                     our spacecraft's orbit propagator and the forces on it
##   4. Attitude                  pointing, sensors, controller, reaction wheels, disturbances
##   5. Visibility and lighting   Earth blocking the view, Earth's shadow, angles
##   6. Saving results            the files written to Outputs/
##   7. Run                       the main time loop
## =================================================================================================
import csv
import json
import os
from dataclasses import dataclass

import numpy as np
import orekit_jpype as orekit

## =================================================================================================
## 1. SETUP
## =================================================================================================

# File locations, relative to this file so it runs from any folder. This file and orekit-data.zip
# are in Propagation/; Inputs/ and Outputs/ are one level up, in the repository folder.
PROPAGATION_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(PROPAGATION_DIR)
INPUT_DIR = os.path.join(REPO_DIR, 'Inputs')
OUTPUT_DIR = os.path.join(REPO_DIR, 'Outputs')

# Start Java and load Orekit's data files. Orekit can only be imported after this.
orekit.initVM()
print('Virtual Machine Initialized...')
from orekit_jpype.pyhelpers import setup_orekit_data
setup_orekit_data(filenames=os.path.join(PROPAGATION_DIR, 'orekit-data.zip'))
print('Orekit data identified...')

from jpype import JArray, JDouble, JImplements, JOverride
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.orekit.forces.drag import DragForce, IsotropicDrag
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel, ThirdBodyAttraction
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient, SolarRadiationPressure
from org.orekit.frames import FramesFactory
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngleType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.propagation.integration import AdditionalDerivativesProvider, CombinedDerivatives
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import AngularCoordinates, Constants, IERSConventions, PVCoordinates

# Frames: gcrf is inertial (every output is in it), itrf rotates with the Earth
gcrf = FramesFactory.getGCRF()
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)  # True: no tidal corrections
utc = TimeScalesFactory.getUTC()

# The Earth (WGS84), Sun and Moon
EARTH_RADIUS_M = Constants.WGS84_EARTH_EQUATORIAL_RADIUS
EARTH_MU = Constants.WGS84_EARTH_MU
earth = OneAxisEllipsoid(EARTH_RADIUS_M, Constants.WGS84_EARTH_FLATTENING, itrf)
sun = CelestialBodyFactory.getSun()
moon = CelestialBodyFactory.getMoon()

## =================================================================================================
## 2. INPUTS
## =================================================================================================

FORCE_MODELS = ('earth_gravity', 'third_body_sun', 'third_body_moon', 'srp', 'drag')
MAX_TLE_AGE_DAYS = 7.0  # warn when a TLE is older than this at the scenario start

@dataclass
class Scenario:
    name: str                 # results go to Outputs/<name>/
    epoch: AbsoluteDate
    duration_days: float
    time_step_s: float
    rso_step_s: float         # how often the RSOs are updated and saved
    rso_file: str             # RSO.json (test set) or GEOTLE.json (full catalog)

@dataclass
class Satellite:
    orbit: KeplerianOrbit     # starting orbit (GCRF)
    true_anomaly_rad: float
    physical: dict            # the physical_properties, propagator and attitude blocks
    propagator: dict
    attitude: dict

@dataclass(frozen=True)
class RSO:
    name: str
    tle: TLE
    object_type: str
    tle_age_days: float       # TLE age at the scenario start

# Contents of one JSON file in Inputs/
def read_input(filename):
    with open(os.path.join(INPUT_DIR, filename), 'r', encoding='utf-8') as file:
        return json.load(file)

# Scenario.json
def read_scenario():
    data = read_input('Scenario.json')
    scenario = Scenario(
        name=data['name'],
        epoch=AbsoluteDate(data['epoch'], utc),
        duration_days=data['duration_days'],
        time_step_s=data['time_step_s'],
        rso_step_s=data['rso_step_s'],
        rso_file=data.get('rso_file', 'RSO.json'),
    )
    steps = scenario.rso_step_s / scenario.time_step_s
    if steps < 1 or abs(steps - round(steps)) > 1e-9:
        raise ValueError('rso_step_s in Inputs/Scenario.json must be a whole multiple '
                         'of time_step_s')
    return scenario

# Satellite.json (the elements are at the scenario start, unless it gives its own "epoch")
def read_satellite(scenario):
    data = read_input('Satellite.json')
    elements = data['orbit']['elements']
    epoch = AbsoluteDate(data['epoch'], utc) if 'epoch' in data else scenario.epoch

    if ('true_anomaly_deg' in elements) == ('start_longitude_deg' in elements):
        raise ValueError('Inputs/Satellite.json orbit elements need exactly one of '
                         'true_anomaly_deg or start_longitude_deg')
    if 'start_longitude_deg' in elements:
        true_anomaly = true_anomaly_at_longitude(elements, epoch, elements['start_longitude_deg'])
    else:
        true_anomaly = np.deg2rad(elements['true_anomaly_deg'])

    unknown_models = [m for m in data['propagator']['force_models'] if m not in FORCE_MODELS]
    if unknown_models:
        raise ValueError(f'unknown force_models {unknown_models} in Inputs/Satellite.json '
                         f'(valid: {", ".join(FORCE_MODELS)})')

    orbit = keplerian_orbit(elements, epoch, true_anomaly)
    return Satellite(orbit, true_anomaly, data['physical_properties'], data['propagator'],
                     data['attitude'])

# The RSO file. Also returns the index of the pointing target.
def read_rsos(scenario, target_name):
    rsos = []
    for entry in read_input(scenario.rso_file)['rsos']:
        line1, line2 = entry['tle_line1'], entry['tle_line2']
        # Orekit silently accepts a bad checksum, so check the format first
        try:
            tle_ok = TLE.isFormatOK(line1, line2)
        except Exception:  # thrown for checksum errors
            tle_ok = False
        if not tle_ok:
            raise ValueError(f"RSO '{entry['name']}' in Inputs/{scenario.rso_file} "
                             "has a malformed TLE")
        tle = TLE(line1, line2)
        age_days = scenario.epoch.durationFrom(tle.getDate()) / 86400.0
        rsos.append(RSO(entry['name'], tle, entry.get('object_type', ''), age_days))

    names = [rso.name for rso in rsos]
    if len(set(names)) != len(names):
        raise ValueError(f'RSO names in Inputs/{scenario.rso_file} must be unique '
                         '(they label the outputs)')

    ages = np.abs([rso.tle_age_days for rso in rsos])
    stale = ages > MAX_TLE_AGE_DAYS
    if stale.any():
        print(f'WARNING: {stale.sum()} of {len(rsos)} TLEs are more than {MAX_TLE_AGE_DAYS:g} days '
              f'from the scenario epoch (worst {ages.max():.1f} days). '
              f'Refresh them for accurate positions.')

    if target_name not in names:
        raise ValueError(f"attitude target_rso '{target_name}' is not in "
                         f"Inputs/{scenario.rso_file}")
    return rsos, names.index(target_name)

# Keplerian orbit (GCRF) from the elements in Satellite.json, at a given true anomaly
def keplerian_orbit(elements, epoch, true_anomaly_rad):
    a = elements['semi_major_axis_km'] * 1000.0
    e = elements['eccentricity']
    i = np.deg2rad(elements['inclination_deg'])
    arg_perigee = np.deg2rad(elements['arg_of_perigee_deg'])
    raan = np.deg2rad(elements['raan_deg'])
    return KeplerianOrbit(a, e, i, arg_perigee, raan, true_anomaly_rad, PositionAngleType.TRUE,
                          gcrf, epoch, EARTH_MU)

# True anomaly that puts the spacecraft over a chosen longitude at the epoch (the Earth turns
# under the orbit, so this changes with the date). Adjusts it until the longitude matches.
def true_anomaly_at_longitude(elements, epoch, longitude_deg):
    target = np.deg2rad(longitude_deg)
    true_anomaly = target
    for _ in range(50):
        orbit = keplerian_orbit(elements, epoch, true_anomaly)
        position = orbit.getPVCoordinates(itrf).getPosition()
        longitude = np.arctan2(position.getY(), position.getX())
        error = (longitude - target + np.pi) % (2 * np.pi) - np.pi  # wrapped to [-pi, pi)
        true_anomaly -= error
        if abs(error) < 1e-12:
            break
    return true_anomaly

## =================================================================================================
## 3. ORBIT
## =================================================================================================

# Our spacecraft's propagator and starting state: Orekit's numerical propagator (Cowell's method,
# Dormand-Prince 8(5,3)) with the force models from Satellite.json, plus the attitude model
def build_propagator(scenario, sat, attitude_dynamics, srp):
    # Steps of at most one time step, 1 m position tolerance, starting each call at a full step
    step = scenario.time_step_s
    tolerances = NumericalPropagator.tolerances(1.0, sat.orbit, sat.orbit.getType())
    integrator = DormandPrince853Integrator(1e-3 * step, step, tolerances[0], tolerances[1])
    integrator.setInitialStepSize(step)

    # Starting state: orbit, mass and attitude [q, w, h] (wheels start with no momentum)
    q = np.array(sat.attitude['initial_quaternion'])
    w = np.deg2rad(sat.attitude['initial_angular_velocity_deg_s'])
    attitude = np.concatenate([q / np.linalg.norm(q), w, np.zeros(3)])
    initial_state = SpacecraftState(sat.orbit, sat.physical['mass_kg'])
    initial_state = initial_state.addAdditionalData(ATTITUDE_STATE, JArray(JDouble)(attitude))

    propagator = NumericalPropagator(integrator)
    propagator.setOrbitType(OrbitType.CARTESIAN)
    propagator.setInitialState(initial_state)
    propagator.addAdditionalDerivativesProvider(attitude_dynamics)

    # Force models (Earth's central gravity is always on)
    models = sat.propagator['force_models']
    if 'earth_gravity' in models:
        field = GravityFieldFactory.getNormalizedProvider(sat.propagator['gravity_degree'],
                                                          sat.propagator['gravity_order'])
        propagator.addForceModel(HolmesFeatherstoneAttractionModel(itrf, field))
    if 'third_body_sun' in models:
        propagator.addForceModel(ThirdBodyAttraction(sun))
    if 'third_body_moon' in models:
        propagator.addForceModel(ThirdBodyAttraction(moon))
    if 'srp' in models:
        propagator.addForceModel(srp)
    if 'drag' in models:
        weather = CssiSpaceWeatherData(CssiSpaceWeatherData.DEFAULT_SUPPORTED_NAMES)
        atmosphere = NRLMSISE00(weather, sun, earth)
        shape = IsotropicDrag(sat.physical['cross_sectional_area_m2'],
                              sat.physical['drag_coefficient'])
        propagator.addForceModel(DragForce(atmosphere, shape))

    return propagator, initial_state

## =================================================================================================
## 4. ATTITUDE
## =================================================================================================

# The attitude is 10 extra numbers Orekit integrates with the orbit: q (quaternion GCRF -> body,
# scalar first), w (body rate [rad/s]) and h (wheel momentum [N m s]). Body +Z (the camera) points
# at the target RSO; body +Y toward the orbit normal (or the zenith when the target is near it).
ATTITUDE_STATE = 'attitude'
BODY_Z = PVCoordinates(Vector3D.PLUS_K, Vector3D.ZERO, Vector3D.ZERO)
BODY_Y = PVCoordinates(Vector3D.PLUS_J, Vector3D.ZERO, Vector3D.ZERO)
ORBIT_NORMAL = 0                     # roll references
ZENITH = 1
ROLL_SWITCH_DEG = (10.0, 15.0)       # to zenith within 10 deg of the orbit normal, back past 15
SRP_PRESSURE_1AU = 4.56e-6           # sunlight pressure at 1 AU [N/m^2]
AU_M = 149597870000.0

# Quaternion product: quat_to_dcm(quat_multiply(p, q)) == quat_to_dcm(p) @ quat_to_dcm(q)
def quat_multiply(p, q):
    scalar = p[0] * q[0] - p[1:] @ q[1:]
    vector = p[0] * q[1:] + q[0] * p[1:] - np.cross(p[1:], q[1:])
    return np.concatenate([[scalar], vector])

# The opposite rotation
def quat_conjugate(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

# Rotation matrix taking vectors from the reference frame into the body frame
def quat_to_dcm(q):
    v = q[1:]
    v_cross = np.array([[0.0, -v[2], v[1]],
                        [v[2], 0.0, -v[0]],
                        [-v[1], v[0], 0.0]])
    return (q[0] ** 2 - v @ v) * np.eye(3) + 2.0 * np.outer(v, v) - 2.0 * q[0] * v_cross

# The attitude model. Orekit calls combinedDerivatives() every integration step; the main loop
# calls update() once per time step.
@JImplements(AdditionalDerivativesProvider)
class AttitudeDynamics:

    def __init__(self, settings, target_propagator, srp_model, physical):
        self.inertia = np.array(settings['inertia_kg_m2'])
        self.inertia_inv = np.linalg.inv(self.inertia)
        self.target_propagator = target_propagator   # its own TLE propagator for the target
        self.roll_reference = ORBIT_NORMAL

        # PD gains from the natural frequency and damping ratio (quaternion error = half angle)
        controller = settings['controller']
        wn = controller['natural_frequency_rad_s']
        self.kp = 2.0 * np.diag(self.inertia) * wn ** 2
        self.kd = 2.0 * controller['damping_ratio'] * wn * np.diag(self.inertia)

        self.max_torque = settings['reaction_wheels']['max_torque_Nm']
        self.max_momentum = settings['reaction_wheels']['max_momentum_Nms']

        # Star tracker and gyro errors, redrawn every time step (seeded, so runs repeat)
        knowledge = settings['knowledge_error']
        self.rng = np.random.default_rng(knowledge['seed'])
        self.attitude_sigma = np.deg2rad(knowledge['attitude_arcsec'] / 3600.0)
        self.rate_sigma = np.deg2rad(knowledge['rate_deg_s'])
        self.attitude_error = np.array([1.0, 0.0, 0.0, 0.0])   # no error until the first update
        self.rate_error = np.zeros(3)

        disturbances = settings['disturbance_torques']
        self.gravity_gradient_on = disturbances['gravity_gradient']
        self.srp_torque_on = disturbances['srp']
        self.center_of_pressure = np.array(disturbances['center_of_pressure_offset_m'])
        self.srp_model = srp_model   # Orekit's sunlight-pressure model, used for Earth's shadow
        self.srp_coefficient = physical['srp_coefficient']
        self.srp_area_m2 = physical['srp_area_m2']

    # Orekit interface (JPype can't call Java default methods, so init and yields are written out)
    @JOverride
    def getName(self):
        return ATTITUDE_STATE

    @JOverride
    def getDimension(self):
        return 10

    @JOverride
    def init(self, initial_state, target):
        pass

    @JOverride
    def yields(self, state):
        return False

    # How fast [q, w, h] are changing (Orekit integrates this)
    @JOverride
    def combinedDerivatives(self, state):
        q, w, wheel_momentum = self.attitude(state)
        q_cmd, w_cmd, w_dot_cmd = self.commanded(state)
        torque = self.wheel_torque(q, w, wheel_momentum, q_cmd, w_cmd, w_dot_cmd)

        # Quaternion kinematics, and Euler's equations with the wheels' momentum
        q_dot = 0.5 * np.concatenate([[-q[1:] @ w], q[0] * w + np.cross(q[1:], w)])
        disturbance = self.disturbance_torque(state, q)
        gyroscopic = np.cross(w, self.inertia @ w + wheel_momentum)
        w_dot = self.inertia_inv @ (torque + disturbance - gyroscopic)
        h_dot = -torque   # the wheels take up the opposite momentum

        derivatives = np.concatenate([q_dot, w_dot, h_dot])
        return CombinedDerivatives(JArray(JDouble)(derivatives), None)

    # Once per time step: new sensor errors, and pick the roll reference
    def update(self, state):
        # Star tracker error (a small random rotation) and gyro error
        half_angles = 0.5 * self.rng.normal(0.0, self.attitude_sigma, 3)
        scalar = np.sqrt(1.0 - half_angles @ half_angles)
        self.attitude_error = np.concatenate([[scalar], half_angles])
        self.rate_error = self.rng.normal(0.0, self.rate_sigma, 3)

        # Roll from the orbit normal, unless the target is too close to that direction
        sat_pv = state.getPVCoordinates(gcrf)
        target = self.target_propagator.getPVCoordinates(state.getDate(), gcrf).getPosition()
        los = np.array(target.subtract(sat_pv.getPosition()).toArray())
        normal = np.array(sat_pv.getMomentum().toArray())
        cross = np.linalg.norm(np.cross(los, normal))
        sin_angle = cross / (np.linalg.norm(los) * np.linalg.norm(normal))

        to_zenith = np.sin(np.deg2rad(ROLL_SWITCH_DEG[0]))
        back_to_normal = np.sin(np.deg2rad(ROLL_SWITCH_DEG[1]))
        if self.roll_reference == ORBIT_NORMAL and sin_angle < to_zenith:
            self.roll_reference = ZENITH
        elif self.roll_reference == ZENITH and sin_angle > back_to_normal:
            self.roll_reference = ORBIT_NORMAL

    # The true attitude: (unit quaternion, body rate, wheel momentum)
    def attitude(self, state):
        x = np.array(state.getAdditionalState(ATTITUDE_STATE))
        q = x[:4] / np.linalg.norm(x[:4])
        return q, x[4:7], x[7:10]

    # The commanded attitude: (quaternion, rate, angular acceleration) from Orekit's
    # AngularCoordinates; the rate and acceleration feed the controller's feedforward
    def commanded(self, state):
        sat_pv = state.getPVCoordinates(gcrf)
        target_pv = self.target_propagator.getPVCoordinates(state.getDate(), gcrf)
        los = PVCoordinates(sat_pv, target_pv).normalize()    # unit vector toward the target

        if self.roll_reference == ORBIT_NORMAL:
            reference = PVCoordinates(sat_pv.getMomentum(), Vector3D.ZERO, Vector3D.ZERO)  # r x v
        else:
            reference = sat_pv                                 # zenith: the position vector

        # Roll reference made perpendicular to the line of sight: (los x ref) x los
        constraint = PVCoordinates.crossProduct(PVCoordinates.crossProduct(los, reference), los)
        constraint = constraint.normalize()

        # Rotation taking the line of sight to body +Z and the constraint to body +Y
        commanded = AngularCoordinates(los, constraint, BODY_Z, BODY_Y, 1.0e-9)
        rotation = commanded.getRotation()
        q_cmd = np.array([rotation.getQ0(), rotation.getQ1(), rotation.getQ2(), rotation.getQ3()])
        w_cmd = np.array(commanded.getRotationRate().toArray())
        w_dot_cmd = np.array(commanded.getRotationAcceleration().toArray())
        return q_cmd, w_cmd, w_dot_cmd

    # Wheel torque on the body [N m]. The controller only sees the measured attitude and rate.
    def wheel_torque(self, q, w, wheel_momentum, q_cmd, w_cmd, w_dot_cmd):
        q_meas = quat_multiply(self.attitude_error, q)   # what the star tracker reports
        w_meas = w + self.rate_error                     # what the gyro reports

        # Attitude error (the short way around) and rate error
        q_err = quat_multiply(q_meas, quat_conjugate(q_cmd))
        if q_err[0] < 0.0:
            q_err = -q_err
        to_body = quat_to_dcm(q_err)
        w_ref = to_body @ w_cmd                  # commanded rate, in body axes
        w_err = w_meas - w_ref

        # Feedforward (follows the commanded motion) plus PD feedback on the error
        feedforward = (self.inertia @ (to_body @ w_dot_cmd - np.cross(w_err, w_ref))
                       + np.cross(w_meas, self.inertia @ w_meas + wheel_momentum))
        proportional = self.kp * q_err[1:]
        derivative = self.kd * w_err
        torque = feedforward - proportional - derivative

        # Wheel torque limit (scaled down, keeping the direction)
        peak = np.max(np.abs(torque))
        if peak > self.max_torque:
            torque *= self.max_torque / peak

        # A wheel at its momentum limit can't spin up any further (no momentum dumping)
        for axis in range(3):
            wheel_full = abs(wheel_momentum[axis]) >= self.max_momentum
            pushing_further = torque[axis] * wheel_momentum[axis] < 0.0
            if wheel_full and pushing_further:
                torque[axis] = 0.0

        return torque

    # Gravity-gradient and sunlight-pressure torques [N m] in body axes, for the true attitude q
    def disturbance_torque(self, state, q):
        to_body = quat_to_dcm(q)
        position = np.array(state.getPosition(gcrf).toArray())
        torque = np.zeros(3)

        # Gravity gradient: 3 mu / r^5 * (r x J r), with r the position in body axes
        if self.gravity_gradient_on:
            r_body = to_body @ position
            scale = 3.0 * EARTH_MU / np.linalg.norm(r_body) ** 5
            torque += scale * np.cross(r_body, self.inertia @ r_body)

        # Sunlight pushes away from the Sun at the center of pressure (lighting: 0 shadow, 1 sun)
        lighting = self.srp_model.getLightingRatio(state) if self.srp_torque_on else 0.0
        if lighting > 0.0:
            to_sun = np.array(sun.getPosition(state.getDate(), gcrf).toArray()) - position
            distance = np.linalg.norm(to_sun)
            pressure = lighting * SRP_PRESSURE_1AU * (AU_M / distance) ** 2
            force_body = (-pressure * self.srp_coefficient * self.srp_area_m2
                          * (to_body @ to_sun) / distance)
            torque += np.cross(self.center_of_pressure, force_body)

        return torque

## =================================================================================================
## 5. VISIBILITY AND LIGHTING
## =================================================================================================

# Positions here are GCRF vectors from the Earth's center [m]
LINE_OF_SIGHT_MARGIN_M = 100e3  # treat the lower atmosphere as blocking the view too

# Angle between vectors [deg] (atan2 stays accurate even for tiny angles)
def angle_deg(a, b):
    sin_part = np.linalg.norm(np.cross(a, b), axis=-1)
    cos_part = np.sum(a * b, axis=-1)
    return np.degrees(np.arctan2(sin_part, cos_part))

# True where the line from the observer to each target passes through the Earth
def earth_blocks_view(observer, targets):
    to_targets = targets - observer
    # Point on each line closest to the Earth's center
    s = np.clip(-(to_targets @ observer) / np.sum(to_targets ** 2, axis=1), 0.0, 1.0)
    closest = np.linalg.norm(observer + s[:, None] * to_targets, axis=1)
    return closest < EARTH_RADIUS_M + LINE_OF_SIGHT_MARGIN_M

# True where a position is in the Earth's shadow (cylinder behind the Earth, one Earth radius wide)
def in_earth_shadow(positions, sun_position):
    sun_direction = sun_position / np.linalg.norm(sun_position)
    along = positions @ sun_direction
    off_axis = np.linalg.norm(positions - along[:, None] * sun_direction, axis=1)
    return (along < 0.0) & (off_axis < EARTH_RADIUS_M)

## =================================================================================================
## 6. SAVING RESULTS
## =================================================================================================

# Output columns (the README describes each one)
POSITION_VELOCITY = ['x_m', 'y_m', 'z_m', 'vx_mps', 'vy_mps', 'vz_mps']
SATELLITE_COLUMNS = ['time_utc', 't_s'] + POSITION_VELOCITY
RSO_COLUMNS = ['time_utc', 't_s', 'name'] + POSITION_VELOCITY
ATTITUDE_COLUMNS = [
    'time_utc', 't_s',
    'q0', 'q1', 'q2', 'q3', 'wx_rad_s', 'wy_rad_s', 'wz_rad_s',
    'q0_cmd', 'q1_cmd', 'q2_cmd', 'q3_cmd',
    'torque_x_Nm', 'torque_y_Nm', 'torque_z_Nm',
    'wheel_hx_Nms', 'wheel_hy_Nms', 'wheel_hz_Nms',
    'disturbance_torque_Nm', 'pointing_error_deg', 'roll_reference',
    'target_range_km', 'target_visible', 'target_sunlit',
    'phase_angle_deg', 'sun_boresight_angle_deg', 'spacecraft_lighting',
]
SUMMARY_COLUMNS = ['name', 'object_type', 'tle_age_days', 'min_range_km', 'min_range_t_s',
                   'visible_pct', 'visible_and_sunlit_pct']

# Opens a CSV file in the output folder and writes its header row
def open_csv(folder, filename, header):
    file = open(os.path.join(folder, filename), 'w', newline='', encoding='utf-8')
    writer = csv.writer(file)
    writer.writerow(header)
    return file, writer

# Position and velocity of an Orekit PVCoordinates, as a list of 6 numbers
def position_velocity(pv):
    return [*pv.getPosition().toArray(), *pv.getVelocity().toArray()]

# initial_state.json: the exact starting state (the STK verification starts from it)
def write_initial_state(folder, sat):
    pv = sat.orbit.getPVCoordinates(gcrf)
    initial_state = {
        'epoch_utc': sat.orbit.getDate().toString(utc),
        'frame': 'GCRF',
        'position_m': list(pv.getPosition().toArray()),
        'velocity_m_s': list(pv.getVelocity().toArray()),
        'true_anomaly_deg': float(np.degrees(sat.true_anomaly_rad) % 360.0),
        'mass_kg': sat.physical['mass_kg'],
    }
    with open(os.path.join(folder, 'initial_state.json'), 'w', encoding='utf-8') as file:
        json.dump(initial_state, file, indent=2)

# rso_summary.csv: each RSO's closest approach and how often it was visible and sunlit
# (sampled every rso_step_s)
class RsoSummary:

    def __init__(self, rsos):
        self.rsos = rsos
        self.min_range_km = np.full(len(rsos), np.inf)
        self.min_range_t_s = np.zeros(len(rsos))
        self.visible_count = np.zeros(len(rsos))
        self.visible_sunlit_count = np.zeros(len(rsos))
        self.samples = 0

    # One sample: every RSO's range [km], visibility and sunlight at time t
    def add(self, t, ranges_km, visible, sunlit):
        closer = ranges_km < self.min_range_km
        self.min_range_km[closer] = ranges_km[closer]
        self.min_range_t_s[closer] = t
        self.visible_count += visible
        self.visible_sunlit_count += visible & sunlit
        self.samples += 1

    # Writes the file (closest first) and prints the five closest approaches
    def write(self, folder):
        order = np.argsort(self.min_range_km)
        ages = np.array([rso.tle_age_days for rso in self.rsos])

        file, writer = open_csv(folder, 'rso_summary.csv', SUMMARY_COLUMNS)
        for i in order:
            writer.writerow([
                self.rsos[i].name,
                self.rsos[i].object_type,
                round(ages[i], 2),
                round(self.min_range_km[i], 1),
                self.min_range_t_s[i],
                round(100.0 * self.visible_count[i] / self.samples, 1),
                round(100.0 * self.visible_sunlit_count[i] / self.samples, 1),
            ])
        file.close()

        print('Closest approaches:')
        for i in order[:5]:
            hours = self.min_range_t_s[i] / 3600.0
            print(f'  {self.rsos[i].name:32s} {self.min_range_km[i]:9.1f} km at t = {hours:6.2f} h')

## =================================================================================================
## 7. RUN
## =================================================================================================

def run():
    scenario = read_scenario()
    sat = read_satellite(scenario)
    rsos, target = read_rsos(scenario, sat.attitude['target_rso'])

    # The RSOs: Orekit's SGP4/SDP4 from their TLEs
    rso_propagators = [TLEPropagator.selectExtrapolator(rso.tle) for rso in rsos]

    # Our spacecraft: sunlight pressure (cannonball), the attitude model (with its own copy of the
    # target's propagator) and the orbit propagator
    surface = IsotropicRadiationSingleCoefficient(sat.physical['srp_area_m2'],
                                                  sat.physical['srp_coefficient'])
    srp = SolarRadiationPressure(sun, earth, surface)
    target_propagator = TLEPropagator.selectExtrapolator(rsos[target].tle)
    dynamics = AttitudeDynamics(sat.attitude, target_propagator, srp, sat.physical)
    propagator, initial_state = build_propagator(scenario, sat, dynamics, srp)

    # Output files
    output_dir = os.path.join(OUTPUT_DIR, scenario.name)
    os.makedirs(output_dir, exist_ok=True)
    write_initial_state(output_dir, sat)
    satellite_file, satellite_csv = open_csv(output_dir, 'satellite.csv', SATELLITE_COLUMNS)
    rsos_file, rsos_csv = open_csv(output_dir, 'rsos.csv', RSO_COLUMNS)
    attitude_file, attitude_csv = open_csv(output_dir, 'attitude.csv', ATTITUDE_COLUMNS)
    summary = RsoSummary(rsos)

    n_steps = int(round(scenario.duration_days * 86400.0 / scenario.time_step_s))
    rso_interval = int(round(scenario.rso_step_s / scenario.time_step_s))
    progress_interval = max(1, int(round(3600.0 / scenario.time_step_s)))  # print hourly
    dynamics.update(initial_state)

    for step in range(n_steps + 1):
        t = step * scenario.time_step_s
        date = scenario.epoch.shiftedBy(t)
        time_utc = date.toString(utc)

        # 1. Our spacecraft (orbit + attitude)
        state = propagator.propagate(date)
        sat_pv = state.getPVCoordinates(gcrf)
        satellite_csv.writerow([time_utc, t] + position_velocity(sat_pv))
        sat_position = np.array(sat_pv.getPosition().toArray())
        sun_position = np.array(sun.getPosition(date, gcrf).toArray())

        # 2. The RSOs: all of them every rso_step_s (saved and summarized), the target every step
        if step % rso_interval == 0:
            rso_pvs = [p.getPVCoordinates(date, gcrf) for p in rso_propagators]
            for rso, pv in zip(rsos, rso_pvs):
                rsos_csv.writerow([time_utc, t, rso.name] + position_velocity(pv))
            positions = np.array([pv.getPosition().toArray() for pv in rso_pvs])
            ranges_km = np.linalg.norm(positions - sat_position, axis=1) / 1e3
            visible = ~earth_blocks_view(sat_position, positions)
            sunlit = ~in_earth_shadow(positions, sun_position)
            summary.add(t, ranges_km, visible, sunlit)
            target_pv = rso_pvs[target]
        else:
            target_pv = rso_propagators[target].getPVCoordinates(date, gcrf)
        target_position = np.array(target_pv.getPosition().toArray())

        # 3. Pointing and lighting
        q, w, wheel_momentum = dynamics.attitude(state)
        q_cmd, w_cmd, w_dot_cmd = dynamics.commanded(state)
        torque = dynamics.wheel_torque(q, w, wheel_momentum, q_cmd, w_cmd, w_dot_cmd)
        disturbance = dynamics.disturbance_torque(state, q)
        los = target_position - sat_position                      # line of sight to the target
        boresight = quat_to_dcm(q).T @ np.array([0.0, 0.0, 1.0])   # body +Z, in GCRF
        pointing_error_deg = angle_deg(boresight, los)
        range_km = np.linalg.norm(los) / 1e3
        target_visible = not earth_blocks_view(sat_position, target_position[None, :])[0]
        target_sunlit = not in_earth_shadow(target_position[None, :], sun_position)[0]
        phase_angle_deg = angle_deg(sun_position - target_position, sat_position - target_position)
        sun_angle_deg = angle_deg(los, sun_position - sat_position)
        lighting = srp.getLightingRatio(state)

        # 4. Save this time step, then new sensor errors and roll reference for the next one
        attitude_csv.writerow([
            time_utc, t, *q, *w, *q_cmd, *torque, *wheel_momentum,
            np.linalg.norm(disturbance), pointing_error_deg, dynamics.roll_reference,
            range_km, int(target_visible), int(target_sunlit),
            phase_angle_deg, sun_angle_deg, lighting,
        ])
        dynamics.update(state)

        if step % progress_interval == 0:
            print(f't = {t / 3600.0:6.2f} h   pointing error = {pointing_error_deg:.6f} deg   '
                  f'target range = {range_km:7.0f} km   '
                  f'wheel momentum = {np.abs(wheel_momentum).max():.4f} N m s')

    # Close the files and write the RSO summary
    for file in (satellite_file, rsos_file, attitude_file):
        file.close()
    summary.write(output_dir)
    print(f'Outputs written to {output_dir}')

if __name__ == '__main__':
    run()
