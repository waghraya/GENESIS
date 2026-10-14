## =================================================================================================
## GENESIS STK VERIFICATION
## Checks the simulation against STK, an independent tool, using the same inputs. It builds the
## matching scenario in STK, exports STK's results, and compares them with the simulation's.
##
## Run the simulation first (python run.py), then on a Windows machine with STK 12 or 13:
##     python Verification/stk_verification.py
## Options:
##     --gravity-file "C:/.../EGM2008.grv"   Earth gravity file for HPOP, if STK's isn't found
##     --rsos NAME NAME ...                  RSOs to build (default: the target + the 4 closest)
##     --compare-only                        skip STK and just re-compare its saved results
##
## Sections:
##   1. Setup              file locations, pass/fail limits, input files, STK names and dates
##   2. Talking to STK     sends STK's Connect commands and logs each one
##   3. Building in STK    the matching STK scenario, propagated and exported
##   4. Reading results    STK's exported files and the simulation's CSVs
##   5. Comparing          positions and attitude, PASS/FAIL
##   6. Run
## =================================================================================================
import argparse
import csv
import glob
import json
import os
import re
import sys
from datetime import datetime, timedelta

import numpy as np

## =================================================================================================
## 1. SETUP
## =================================================================================================

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(REPO_DIR, 'Inputs')
RESULTS_DIR = os.path.join(REPO_DIR, 'Verification', 'STK_Results')   # STK's exports go here
STK_SCENARIO = 'GENESIS_Verification'
MAX_STK_RSOS = 5   # the target plus the closest others (a whole catalog would make STK very slow)

# Pass/fail limits
SAT_POSITION_TOL_M = 100.0    # same force models in both tools; kilometers means one is off
RSO_POSITION_TOL_M = 1000.0   # same TLEs, but different SGP4 implementations and frame conversions
ATTITUDE_TOL_DEG = 0.01
SENSOR_SIGMAS = 6.0           # allowance for the simulated sensor errors in the true attitude
SETTLE_TIME_S = 600.0         # judge the true attitude only after the initial slew
Y_AXIS = np.array([0.0, 1.0, 0.0])
Z_AXIS = np.array([0.0, 0.0, 1.0])

MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Contents of one JSON file in Inputs/
def read_input(filename):
    with open(os.path.join(INPUT_DIR, filename), encoding='utf-8') as file:
        return json.load(file)

# The RSO entries the scenario uses (its rso_file, RSO.json if not set)
def read_rsos(scenario):
    return read_input(scenario.get('rso_file', 'RSO.json'))['rsos']

# STK object names only allow letters, digits and underscores
def stk_name(name):
    return re.sub(r'[^A-Za-z0-9_]', '_', name)

# Our spacecraft's STK name (from Satellite.json), also the name of its .e/.a result files
SAT_NAME = stk_name(read_input('Satellite.json')['satellite']['name'])
SAT_PATH = f'*/Satellite/{SAT_NAME}'

# ISO time ('2026-10-14T00:00:00Z', with or without the Z) -> datetime in UTC
def parse_iso(text):
    return datetime.fromisoformat(text.replace('Z', '+00:00')).replace(tzinfo=None)

# datetime -> STK's date format, e.g. '1 Jan 2026 00:00:00.000'
def stk_date(dt):
    return f'{dt.day} {MONTHS[dt.month - 1]} {dt.year} {dt:%H:%M:%S}.{dt.microsecond // 1000:03d}'

# STK's date format ('1 Jan 2026 00:00:00.000000') -> datetime
def parse_stk_date(text):
    day, month, year, clock = text.split()[:4]
    hours, minutes, seconds = clock.split(':')
    whole, _, fraction = seconds.partition('.')
    return datetime(int(year), MONTHS.index(month[:3]) + 1, int(day), int(hours), int(minutes),
                    int(whole), int((fraction + '000000')[:6]))

## =================================================================================================
## 2. TALKING TO STK
## =================================================================================================

# Attaches to (or starts) STK and returns a function that runs one Connect command (STK's text
# command language, which works the same across STK versions)
def connect_to_stk():
    try:
        from ansys.stk.core.stkdesktop import STKDesktop   # STK 12.10+ and 13
        try:
            stk = STKDesktop.attach_to_application()
        except Exception:
            stk = STKDesktop.start_application(visible=True, user_control=True)
        return stk.root.execute_command
    except ImportError:
        from agi.stk12.stkdesktop import STKDesktop        # STK 12.0 - 12.9
        try:
            stk = STKDesktop.AttachToApplication()
        except Exception:
            stk = STKDesktop.StartApplication(visible=True, userControl=True)
        return stk.Root.ExecuteCommand

# Sends Connect commands and logs each one to stk_commands.log, with whether STK accepted it
class Connect:

    def __init__(self, execute, log_path):
        self.execute = execute
        self.log = open(log_path, 'w', encoding='utf-8')

    # STK raises an error when it rejects a command: a required command stops the script, an
    # optional one just returns False
    def run(self, command, required=True):
        print(f'> {command}')
        try:
            self.execute(command)
        except Exception as err:
            self.log.write(f'REJECTED  {command}\n          {err}\n')
            self.log.flush()
            if required:
                raise RuntimeError(f'STK rejected a required command:\n  {command}\n'
                                   f'  {err}') from err
            print('  (rejected, continuing)')
            return False
        self.log.write(f'OK        {command}\n')
        self.log.flush()
        return True

    # Some command options differ between STK versions: runs the first variant STK accepts
    def run_first(self, commands, description):
        for command in commands:
            if self.run(command, required=False):
                return command
        raise RuntimeError(f'STK rejected every variant for {description}; see {self.log.name}')

## =================================================================================================
## 3. BUILDING IN STK
## =================================================================================================

# The pointing target plus the closest others from the simulation's rso_summary.csv
def select_rsos(rsos, target, summary_path, requested):
    by_name = {rso['name']: rso for rso in rsos}
    if requested:
        unknown = [name for name in requested if name not in by_name]
        if unknown:
            raise ValueError(f'--rsos names not in the scenario\'s RSO file: {unknown}')
        names = list(requested)
    else:
        others = [rso['name'] for rso in rsos]
        if os.path.exists(summary_path):
            with open(summary_path, newline='', encoding='utf-8') as file:
                # rso_summary.csv lists the closest first
                others = [row['name'] for row in csv.DictReader(file) if row['name'] in by_name]
        names = [target] + [name for name in others if name != target][:MAX_STK_RSOS - 1]
    if target not in names:
        names.insert(0, target)
    return [by_name[name] for name in names]

# Clears old exports (so different runs never mix) and writes the TLE file STK reads
def prepare_results_dir(rsos):
    os.makedirs(RESULTS_DIR, exist_ok=True)
    for pattern in ('*.e', '*.a'):
        for old_file in glob.glob(os.path.join(RESULTS_DIR, pattern)):
            os.remove(old_file)
    tle_path = os.path.join(RESULTS_DIR, 'rso.tle')
    with open(tle_path, 'w', encoding='ascii') as file:
        for rso in rsos:
            file.write(f"{rso['tle_line1']}\n{rso['tle_line2']}\n")
    return tle_path

# Looks for an Earth gravity file in the usual STK install folders (the newest STK if several)
def find_gravity_file():
    for name in ('EGM2008.grv', 'WGS84_EGM96.grv', 'EGM96.grv'):
        for root in (r'C:\Program Files', r'C:\Program Files (x86)'):
            for depth in ('*', os.path.join('*', '*'), os.path.join('*', '*', '*')):
                matches = glob.glob(os.path.join(root, depth, 'STKData', 'CentralBodies', 'Earth',
                                                 name))
                if matches:
                    return sorted(matches)[-1]
    return None

def create_scenario(connect, start, interval):
    connect.run('Unload / *', required=False)   # close whatever scenario is open
    connect.run(f'New / Scenario {STK_SCENARIO}')
    connect.run(f'SetAnalysisTimePeriod * {interval}')
    connect.run(f'SetEpoch * "{stk_date(start)}"', required=False)
    connect.run('SetUnits / m', required=False)   # already the default, just making sure

# HPOP force models matching the force_models and physical_properties in Satellite.json
def configure_forces(connect, satellite, gravity_file):
    forces = satellite['propagator']['force_models']
    physical = satellite['physical_properties']

    if 'earth_gravity' in forces:
        degree = satellite['propagator']['gravity_degree']
        order = satellite['propagator']['gravity_order']
    else:
        degree, order = 0, 0   # point mass only
    gravity_command = f'HPOP {SAT_PATH} Force Gravity "{gravity_file}" {degree} {order}'
    if not gravity_file or not connect.run(gravity_command, required=False):
        # Mostly harmless at GEO: terms above degree ~4 are below 1e-12 m/s^2 there
        print(f'WARNING: could not set HPOP gravity to {degree}x{order}, '
              f'keeping the STK default gravity field')

    for body, model in (('Sun', 'third_body_sun'), ('Moon', 'third_body_moon')):
        if model in forces:
            # JPLDEFile takes the body's gravity constant from the JPL tables, like Orekit does
            commands = [f'HPOP {SAT_PATH} Force ThirdBodyGravity {body} On {source}'
                        for source in ('JPLDEFile', 'FromCB', 'Internal')]
            connect.run_first(commands, f'{body} gravity')
        else:
            connect.run(f'HPOP {SAT_PATH} Force ThirdBodyGravity {body} Off')

    if 'srp' in forces:
        cr = physical['srp_coefficient']
        area_mass = physical['srp_area_m2'] / physical['mass_kg']
        connect.run_first([
            f'HPOP {SAT_PATH} Force SRP On Model Spherical {cr} {area_mass} ShadowModel DualCone',
            f'HPOP {SAT_PATH} Force SRP On Model Spherical {cr} {area_mass}',
        ], 'solar radiation pressure')
        # The simulation only uses Earth's shadow
        connect.run(f'HPOP {SAT_PATH} Force EclipsingBodies Remove Moon', required=False)
    else:
        connect.run(f'HPOP {SAT_PATH} Force SRP Off')

    if 'drag' in forces:
        # A different density model than the simulation's, but density is ~zero at GEO either way
        drag_area_mass = physical['cross_sectional_area_m2'] / physical['mass_kg']
        connect.run(f'HPOP {SAT_PATH} Drag On {physical["drag_coefficient"]} {drag_area_mass} '
                    f'"Jacchia 1960"', required=False)
    else:
        connect.run(f'HPOP {SAT_PATH} Drag Off', required=False)
    connect.run(f'HPOP {SAT_PATH} Force SolidTides Off', required=False)   # not in the simulation

# Our spacecraft, propagated with HPOP from the simulation's exact starting state (the first
# SetState switches it to HPOP so the forces can be edited, the second re-propagates)
def add_spacecraft(connect, satellite, initial_state, interval, step, gravity_file):
    values = initial_state['position_m'] + initial_state['velocity_m_s']
    state = ' '.join(f'{value:.9f}' for value in values)
    epoch = stk_date(parse_iso(initial_state['epoch_utc']))
    set_state = f'SetState {SAT_PATH} Cartesian HPOP {interval} {step} ICRF "{epoch}" {state}'
    connect.run(f'New / */Satellite {SAT_NAME}')
    connect.run(set_state)
    configure_forces(connect, satellite, gravity_file)
    connect.run(set_state)

# One STK satellite per RSO, propagated with SGP4 from the same TLE
def add_rsos(connect, rsos, interval, step, tle_path):
    for rso in rsos:
        name = stk_name(rso['name'])
        ssc = rso['tle_line1'][2:7].strip()   # catalog number, as written in the TLE
        connect.run(f'New / */Satellite {name}')
        connect.run_first([
            f'SetState */Satellite/{name} SGP4 {interval} {step} {ssc} '
            f'TLESource Automatic Source File "{tle_path}"',
            f'SetState */Satellite/{name} TLE "{rso["tle_line1"]}" "{rso["tle_line2"]}" '
            f'TimePeriod {interval} Routine SGP4',
        ], f'SGP4 propagation of {name}')

# Aligned and constrained attitude: body +Z toward the target, body +Y toward the orbit normal.
# Vector names vary between STK versions, so several are tried.
def set_attitude(connect, target):
    target_name = stk_name(target)
    connect.run(f'VectorTool * Satellite/{SAT_NAME} Create Vector ToTarget Displacement '
                f'"Satellite/{SAT_NAME} Center" "Satellite/{target_name} Center"', required=False)
    align_vectors = [
        'ToTarget',
        f'Satellite/{SAT_NAME} ToTarget',
        f'Satellite/{SAT_NAME} {target_name}',
    ]
    constraint_vectors = ['Orbit_Normal', 'OrbitNormal', 'Orbit_AngMomentum', 'OrbitMomentum']
    connect.run_first([
        f'SetAttitude {SAT_PATH} Profile AlignConstrain Axis 0 0 1 "{align}" '
        f'Axis 0 1 0 "{constraint}"'
        for align in align_vectors for constraint in constraint_vectors
    ], 'the aligned/constrained attitude profile')

# Exports every object's position and velocity (ICRF) and our attitude, at the scenario time step
def export_results(connect, rsos, interval, step):
    for name in [SAT_NAME] + [stk_name(rso['name']) for rso in rsos]:
        path = os.path.join(RESULTS_DIR, name + '.e')
        connect.run(f'ExportDataFile */Satellite/{name} Ephemeris "{path}" Type STK CoordSys ICRF '
                    f'CentralBody Earth TimeSteps {step} TimePeriod {interval}')
    attitude_file = os.path.join(RESULTS_DIR, SAT_NAME + '.a')
    options = f'Details Quaternions TimeSteps {step} TimePeriod {interval}'
    connect.run_first([
        f'ExportDataFile {SAT_PATH} Attitude "{attitude_file}" '
        f'CoordAxes Custom "CentralBody/Earth ICRF" {options}',
        f'ExportDataFile {SAT_PATH} Attitude "{attitude_file}" CoordAxes J2000 {options}',
    ], 'the attitude export')

# Builds the matching scenario in STK, propagates it and exports the results
def build_in_stk(scenario, satellite, sim_dir, gravity_file, requested_rsos):
    with open(os.path.join(sim_dir, 'initial_state.json'), encoding='utf-8') as file:
        initial_state = json.load(file)   # the exact starting state the simulation used
    start = parse_iso(scenario['epoch'])
    stop = start + timedelta(days=scenario['duration_days'])
    step = scenario['time_step_s']
    interval = f'"{stk_date(start)}" "{stk_date(stop)}"'
    target = satellite['attitude']['target_rso']
    summary_path = os.path.join(sim_dir, 'rso_summary.csv')
    rsos = select_rsos(read_rsos(scenario), target, summary_path, requested_rsos)
    tle_path = prepare_results_dir(rsos)

    connect = Connect(connect_to_stk(), os.path.join(RESULTS_DIR, 'stk_commands.log'))
    create_scenario(connect, start, interval)
    gravity_file = gravity_file or find_gravity_file()
    add_spacecraft(connect, satellite, initial_state, interval, step, gravity_file)
    add_rsos(connect, rsos, interval, step, tle_path)
    set_attitude(connect, target)
    export_results(connect, rsos, interval, step)
    print(f'\nSTK results for {SAT_NAME} and {len(rsos)} RSOs written to {RESULTS_DIR}\n')

## =================================================================================================
## 4. READING RESULTS
## =================================================================================================

# Data block and header of an STK .e (positions) or .a (quaternions) file. Times are shifted to
# seconds from the scenario epoch, and distances converted to meters.
def read_stk_file(path, data_keyword, scenario_epoch):
    stk_epoch = None
    header = {}
    rows = []
    in_data = False
    distance_scale = 1.0
    with open(path, encoding='utf-8', errors='replace') as file:
        for line in file:
            tokens = line.split()
            if not tokens:
                continue
            if in_data:
                if tokens[0] == 'END':
                    break
                rows.append([float(value) for value in tokens])
            elif tokens[0] == data_keyword:
                in_data = True
            else:
                header[tokens[0]] = ' '.join(tokens[1:])
                if tokens[0] == 'ScenarioEpoch':
                    stk_epoch = parse_stk_date(header['ScenarioEpoch'])
                elif tokens[0] == 'DistanceUnit' and tokens[1].lower().startswith('kilomet'):
                    distance_scale = 1000.0
    if stk_epoch is None or not rows:
        raise ValueError(f'{path}: could not find ScenarioEpoch and a {data_keyword} block')
    data = np.array(rows)
    data[:, 0] += (stk_epoch - scenario_epoch).total_seconds()
    if data_keyword == 'EphemerisTimePosVel':
        data[:, 1:7] *= distance_scale
    return data, header

# Numeric columns of a simulation CSV (skipping the time_utc text column), and their names
def read_sim_csv(path):
    with open(path, newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        header = next(reader)
        rows = [[float(value) for value in row[1:]] for row in reader]
    return np.array(rows), header[1:]

# Rows of rsos.csv for the given RSOs: {name: array of [t_s, x, y, z, vx, vy, vz]}
def read_sim_rsos(path, names):
    rows = {name: [] for name in names}
    with open(path, newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            if row[2] in rows:
                rows[row[2]].append([float(row[1])] + [float(value) for value in row[3:9]])
    return {name: np.array(data) for name, data in rows.items() if data}

# Indices of the times both series have (to the millisecond)
def match_times(t_a, t_b):
    index_b = {round(t, 3): i for i, t in enumerate(t_b)}
    pairs = [(i, index_b[round(t, 3)]) for i, t in enumerate(t_a) if round(t, 3) in index_b]
    if not pairs:
        raise ValueError('the simulation and STK outputs have no epochs in common')
    return np.array(pairs).T

# Indices into times of the given epochs (to the millisecond)
def select_epochs(times, epochs):
    index = {round(t, 3): i for i, t in enumerate(times)}
    return np.array([index[epoch] for epoch in epochs])

## =================================================================================================
## 5. COMPARING
## =================================================================================================

# Rotation matrices (reference -> body) for scalar-first quaternions, shape (N, 4) -> (N, 3, 3)
def dcm(q):
    q0, v = q[:, 0], q[:, 1:]
    v_cross = np.zeros((len(q), 3, 3))
    v_cross[:, 0, 1], v_cross[:, 0, 2], v_cross[:, 1, 2] = -v[:, 2], v[:, 1], -v[:, 0]
    v_cross -= v_cross.transpose(0, 2, 1)
    return ((q0 ** 2 - np.sum(v * v, axis=1))[:, None, None] * np.eye(3)
            + 2.0 * v[:, :, None] * v[:, None, :] - 2.0 * q0[:, None, None] * v_cross)

# A body axis (e.g. +Z) in the reference frame, for each quaternion
def body_axis_in_reference(q, axis):
    return np.einsum('nji,j->ni', dcm(q), axis)   # transpose(dcm) @ axis

# Angle between matching rows of a and b [deg]
def angle_deg(a, b):
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.sum(a * b, axis=1)))

def unit(vectors):
    return vectors / np.linalg.norm(vectors, axis=1)[:, None]

# Prints one line per check and remembers whether any failed
class Report:

    def __init__(self):
        self.failed = False

    # values: one per epoch; limit: a number, or one limit per epoch
    def check(self, label, values, limit, unit_label, limit_label=None):
        ok = bool(np.all(values <= limit))
        self.failed |= not ok
        status = 'PASS' if ok else 'FAIL'
        rms = np.sqrt(np.mean(values ** 2))
        limit_text = limit_label or f'{limit:g}'
        print(f'  [{status}] {label:<46} max {np.max(values):11.4g} {unit_label:<4}'
              f' rms {rms:11.4g}   final {values[-1]:11.4g}   (limit {limit_text})')

# Position check. sim is an array of [t_s, x, y, z, vx, vy, vz]; returns the common epochs and
# both states there.
def compare_ephemeris(report, label, sim, stk_path, scenario_epoch, limit_m):
    stk, header = read_stk_file(stk_path, 'EphemerisTimePosVel', scenario_epoch)
    i_sim, i_stk = match_times(sim[:, 0], stk[:, 0])
    if header.get('CoordinateSystem', 'ICRF') != 'ICRF':
        print(f'  note: {os.path.basename(stk_path)} is in {header["CoordinateSystem"]}, not ICRF')
    position_error = np.linalg.norm(sim[i_sim, 1:4] - stk[i_stk, 1:4], axis=1)
    velocity_error = np.linalg.norm(sim[i_sim, 4:7] - stk[i_stk, 4:7], axis=1)
    report.check(f'{label} position', position_error, limit_m, 'm')
    print(f'         {label} velocity difference: max {velocity_error.max():.4g} m/s '
          f'over {len(i_sim)} epochs')
    return sim[i_sim, 0], sim[i_sim, 1:7], stk[i_stk, 1:7]

# The three attitude checks. sat_states and target_states are (epochs, sim states, STK states).
def compare_attitude(report, scenario_epoch, attitude_inputs, sim_dir, sat_states, target_states):
    print(f'\nAttitude (target {attitude_inputs["target_rso"]})')
    sim, columns = read_sim_csv(os.path.join(sim_dir, 'attitude.csv'))
    stk_path = os.path.join(RESULTS_DIR, f'{SAT_NAME}.a')
    stk, header = read_stk_file(stk_path, 'AttitudeTimeQuaternions', scenario_epoch)
    if header.get('CoordinateAxes', 'ICRF') != 'ICRF':
        # J2000 and ICRF differ by ~0.02 arcsec
        print(f'  note: STK attitude is relative to {header["CoordinateAxes"]} axes')

    # Epochs where everything exists (the RSOs are only saved every rso_step_s)
    t_sat, sim_sat, stk_sat = sat_states
    t_target, sim_target, stk_target = target_states
    common = sorted({round(t, 3) for t in sim[:, 0]} & {round(t, 3) for t in stk[:, 0]}
                    & {round(t, 3) for t in t_sat} & {round(t, 3) for t in t_target})
    sim = sim[select_epochs(sim[:, 0], common)]
    stk = stk[select_epochs(stk[:, 0], common)]
    sat_sim = sim_sat[select_epochs(t_sat, common)]
    sat_stk = stk_sat[select_epochs(t_sat, common)]
    target_sim = sim_target[select_epochs(t_target, common)]
    target_stk = stk_target[select_epochs(t_target, common)]
    t = sim[:, 0]

    los_sim = unit(target_sim[:, :3] - sat_sim[:, :3])
    los_stk = unit(target_stk[:, :3] - sat_stk[:, :3])
    q_stk = stk[:, [4, 1, 2, 3]]   # STK puts the scalar last
    q_true = sim[:, [columns.index(c) for c in ('q0', 'q1', 'q2', 'q3')]]
    q_cmd = sim[:, [columns.index(c) for c in ('q0_cmd', 'q1_cmd', 'q2_cmd', 'q3_cmd')]]
    orbit_normal_roll = sim[:, columns.index('roll_reference')] == 0   # STK always uses it

    # LoS offset: the angle between the sim's and STK's lines of sight (from position differences).
    # +Y is the orbit normal projected off the line of sight, so it can move (1 + |cot|) times that.
    los_offset = angle_deg(los_sim, los_stk)
    orbit_normal = unit(np.cross(sat_sim[:, :3], sat_sim[:, 3:6]))
    to_normal = np.radians(angle_deg(los_sim, orbit_normal))
    y_offset = los_offset * (1.0 + np.abs(np.cos(to_normal) / np.sin(to_normal)))
    print(f'  {len(t)} epochs; sim and STK lines of sight differ by up to '
          f'{los_offset.max():.4g} deg (from the position differences)')
    with_offset = f'{ATTITUDE_TOL_DEG:g} + LoS offset'

    # 1) Convention check: does STK's own +Z sit on its line of sight? (hint if the conjugate fits)
    stk_z = body_axis_in_reference(q_stk, Z_AXIS)
    stk_boresight = angle_deg(stk_z, los_stk)
    report.check('STK +Z vs STK line of sight (convention)', stk_boresight, ATTITUDE_TOL_DEG, 'deg')
    if stk_boresight.max() > ATTITUDE_TOL_DEG:
        flipped = angle_deg(body_axis_in_reference(q_stk * [1, -1, -1, -1], Z_AXIS), los_stk)
        cause = ('quaternion convention mismatch' if flipped.max() <= ATTITUDE_TOL_DEG
                 else 'geometry mismatch')
        print(f'         with the conjugate quaternion this would be {flipped.max():.4g} deg '
              f'({cause})')

    # 2) Commanded attitude; roll only where the sim also used the orbit normal
    cmd_z = angle_deg(body_axis_in_reference(q_cmd, Z_AXIS), stk_z)
    report.check('Commanded +Z vs STK +Z', cmd_z, ATTITUDE_TOL_DEG + los_offset, 'deg', with_offset)
    if orbit_normal_roll.any():
        cmd_y = angle_deg(body_axis_in_reference(q_cmd[orbit_normal_roll], Y_AXIS),
                          body_axis_in_reference(q_stk[orbit_normal_roll], Y_AXIS))
        report.check('Commanded +Y vs STK +Y (roll)', cmd_y,
                     ATTITUDE_TOL_DEG + y_offset[orbit_normal_roll], 'deg', with_offset)
    if not orbit_normal_roll.all():
        print(f'         roll not compared at {np.sum(~orbit_normal_roll)} epochs '
              f'where the sim used its zenith roll reference (target near the orbit normal)')

    # 3) True attitude. Sensor allowance: the attitude error plus the offset the controller
    #    settles at while chasing a rate error (2 zeta / wn times the rate error)
    knowledge = attitude_inputs['knowledge_error']
    controller = attitude_inputs['controller']
    sensor_deg = SENSOR_SIGMAS * (knowledge['attitude_arcsec'] / 3600.0
                                  + 2.0 * controller['damping_ratio'] * knowledge['rate_deg_s']
                                  / controller['natural_frequency_rad_s'])
    settled = t >= SETTLE_TIME_S
    true_boresight = angle_deg(body_axis_in_reference(q_true[settled], Z_AXIS), los_stk[settled])
    report.check(f'True +Z vs STK line of sight (t >= {SETTLE_TIME_S:g} s)', true_boresight,
                 ATTITUDE_TOL_DEG + los_offset[settled] + sensor_deg, 'deg',
                 f'{with_offset} + sensor error')

    # How long the initial slew took, from the sim's own pointing error
    pointing = sim[:, columns.index('pointing_error_deg')]
    limit = ATTITUDE_TOL_DEG + sensor_deg
    above = np.nonzero(pointing > limit)[0]
    if len(above) and above[-1] == len(t) - 1:
        settle = 'never settles'
    else:
        settle_time = t[above[-1] + 1] if len(above) else t[0]
        settle = f'below {limit:.3g} deg after t = {settle_time:g} s'
    print(f'         simulated slew: {pointing[0]:.2f} deg at t = {t[0]:g} s, {settle}')

# Compares the simulation with STK's results (only the RSOs STK exported); returns 1 on a FAIL
def compare(scenario, attitude_inputs, sim_dir):
    scenario_epoch = parse_iso(scenario['epoch'])
    target_name = attitude_inputs['target_rso']
    report = Report()

    if not os.path.isdir(RESULTS_DIR):
        raise FileNotFoundError(f'{RESULTS_DIR} not found: run without --compare-only first')
    stk_files = {os.path.splitext(f)[0] for f in os.listdir(RESULTS_DIR)
                 if f.endswith('.e') and f != f'{SAT_NAME}.e'}
    compared = [rso['name'] for rso in read_rsos(scenario) if stk_name(rso['name']) in stk_files]
    if target_name not in compared:
        raise ValueError(f'{RESULTS_DIR} has no ephemeris for the pointing target {target_name}')
    sim_rsos = read_sim_rsos(os.path.join(sim_dir, 'rsos.csv'), compared)

    print(f'Simulation: {sim_dir}\nSTK:        {RESULTS_DIR}\n')
    print('Ephemerides (GCRF vs STK ICRF)')
    sim_sat, _ = read_sim_csv(os.path.join(sim_dir, 'satellite.csv'))
    sat_states = compare_ephemeris(report, 'Spacecraft', sim_sat,
                                   os.path.join(RESULTS_DIR, f'{SAT_NAME}.e'),
                                   scenario_epoch, SAT_POSITION_TOL_M)
    rso_states = {}
    for name in compared:
        path = os.path.join(RESULTS_DIR, f'{stk_name(name)}.e')
        rso_states[name] = compare_ephemeris(report, name, sim_rsos[name], path, scenario_epoch,
                                             RSO_POSITION_TOL_M)

    compare_attitude(report, scenario_epoch, attitude_inputs, sim_dir, sat_states,
                     rso_states[target_name])

    print('\nRESULT:', 'FAIL' if report.failed else 'PASS')
    return 1 if report.failed else 0

## =================================================================================================
## 6. RUN
## =================================================================================================

def main():
    parser = argparse.ArgumentParser(description='Checks the simulation against STK.')
    parser.add_argument('--gravity-file', help='Earth .grv file for HPOP (default: search the '
                                               'STK install)')
    parser.add_argument('--rsos', nargs='+', help=f'RSO names to build (default: the target plus '
                                                  f'the closest, {MAX_STK_RSOS} in total)')
    parser.add_argument('--compare-only', action='store_true',
                        help='skip STK and compare the STK results already saved')
    args = parser.parse_args()

    # The same inputs as the simulation, and the results of its last run
    scenario = read_input('Scenario.json')
    satellite = read_input('Satellite.json')
    sim_dir = os.path.join(REPO_DIR, 'Outputs', scenario['name'])
    if not os.path.exists(os.path.join(sim_dir, 'initial_state.json')):
        raise FileNotFoundError(f'no simulation results in {sim_dir}: run python run.py first')

    if not args.compare_only:
        build_in_stk(scenario, satellite, sim_dir, args.gravity_file, args.rsos)
    return compare(scenario, satellite['attitude'], sim_dir)

if __name__ == '__main__':
    sys.exit(main())
