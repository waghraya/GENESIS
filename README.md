# GENESIS

Mission simulator for a small GEO spacecraft that observes nearby resident space objects (RSOs): active satellites, dead satellites, rocket bodies and debris. It propagates the spacecraft's orbit and attitude and the RSOs around it, and checks pointing, visibility and lighting at every step. The repository also holds the tools to load real satellite data from Space-Track and to verify the simulator against STK.

**Contents:** [Quick start](#quick-start) · [Repository layout](#repository-layout) · [How the simulation works](#how-the-simulation-works) · [Inputs](#inputs) · [Outputs](#outputs) · [Updating the GEO catalog](#updating-the-geo-catalog-spacetrack) · [STK verification](#stk-verification) · [OrbitCalculations](#orbitcalculations) · [Status](#status)

## Quick start

1. Install a **Java JDK 11+** and set the `JAVA_HOME` environment variable. The simulator runs [Orekit](https://www.orekit.org/), a Java library, through the `orekit_jpype` Python wrapper.
2. Install the Python packages (Python 3.13): `pip install -r requirements.txt`
3. Run `python run.py`, from any folder. It reads `Inputs/` and writes `Outputs/<scenario name>/`. The 3-day test case takes about 3–4 minutes.

## Repository layout

```
GENESIS/
├── run.py              the start button: runs the simulation
├── requirements.txt    Python packages
├── Propagation/        the simulation
│   ├── propagator.py       the simulation itself (one script)
│   └── orekit-data.zip     data Orekit needs (Earth orientation, ephemerides, gravity field, space weather)
├── Inputs/             what to simulate: scenario, spacecraft, RSOs
├── Outputs/            simulation results (created by run.py, not in git)
├── SpaceTrack/         turns a Space-Track download into the GEO catalog
├── Verification/       checks the simulator against STK (stk_verification.py)
└── OrbitCalculations/  MATLAB mission analysis (SRR)
```

## How the simulation works

**`run.py`** is the start button. It calls the main loop in `Propagation/propagator.py` and prints "Sim finished!" at the end.

**`Propagation/propagator.py`** is the whole simulation, in numbered sections that follow the order of a run. Orekit propagates both our spacecraft and the RSOs:
- **our spacecraft's orbit:** Orekit's numerical propagator with Orekit's force models;
- **the RSOs:** Orekit's SGP4/SDP4 from their TLEs;
- **our spacecraft's attitude:** our own equations (Orekit has no controller, wheels or torques), integrated by Orekit alongside the orbit.

At each time step it moves everything forward and records pointing, visibility and lighting.

| Section | What it does |
|---|---|
| **1. Setup** | Starts Java and loads `orekit-data.zip`. Sets up the reference frames (GCRF, ITRF), UTC, and Earth, Sun and Moon. |
| **2. Inputs** | Reads and checks the three input files:<ul><li>builds the starting orbit, solving for the start longitude if one is given;</li><li>validates the RSO TLEs and warns about stale ones.</li></ul> |
| **3. Orbit** | The spacecraft's numerical propagator: Cowell's method with a Dormand-Prince 8(5,3) integrator, plus the force models listed in `Satellite.json` (Earth gravity field, Sun/Moon gravity, solar radiation pressure, drag) |
| **4. Attitude** | The attitude model:<ul><li>pointing: camera (+Z) on the target, roll from the orbit normal (or the zenith near it);</li><li>star tracker and gyro errors;</li><li>PD controller with feedforward;</li><li>reaction wheels with torque and momentum limits;</li><li>gravity-gradient and solar pressure torques.</li></ul> |
| **5. Visibility and lighting** | Whether Earth blocks the view, whether objects are in Earth's shadow, and angles |
| **6. Saving results** | The output columns, `initial_state.json`, and the per-RSO summary (`rso_summary.csv`) |
| **7. Run** | The main loop that uses all of the above |

**Conventions:**
- **Units:** SI (m, m/s, rad internally; degrees in inputs and outputs where the name says `_deg`). Units are in the field and column names.
- **Frame:** all positions and attitudes are in **GCRF**, the inertial frame (STK calls it ICRF).
- **Time:** UTC.
- **Quaternions:** scalar first, rotation GCRF → body.

## Inputs

Four JSON files in `Inputs/`. JSON can't hold comments, so every field is explained here.

| File | Contents |
|---|---|
| `Scenario.json` | Start time, duration, time steps, and which RSO list to use |
| `Satellite.json` | Orbit, physical properties, force models, attitude/controller/wheel/sensor settings |
| `RSO.json` | A small RSO list (name + TLE): the 3-object test set |
| `GEOTLE.json` | The full GEO catalog from Space-Track (not in git) |

The current case uses the MATLAB disposal script's initial orbit, starting at 114.7°W on Oct 14, 2026, and points at a debris piece. **The physical, inertia, controller, wheel, sensor and center-of-pressure values are placeholders** until the design has real numbers.

### Scenario.json

| Field | Meaning |
|---|---|
| `name` | Output folder: `Outputs/<name>/` |
| `epoch` | Start time, UTC |
| `duration_days` | Length of the run |
| `time_step_s` | Step for the spacecraft orbit, attitude and the pointing target |
| `rso_step_s` | How often every RSO is propagated, saved to `rsos.csv` and counted in `rso_summary.csv`. Must be a whole multiple of `time_step_s`. Raise it (e.g. 600) before running the full catalog. |
| `rso_file` | Which RSO list to use: `RSO.json` (the test set) or `GEOTLE.json` (the full catalog). Defaults to `RSO.json`. |

### Satellite.json

- **`satellite`:** `name` is used as the spacecraft's object name in STK and for its STK result files. `id` and `international_designator` are kept for when the spacecraft has catalog numbers.
- **`epoch`** (optional): the time the orbital elements are given at. If it's left out, the scenario's `epoch` is used, so only `Scenario.json` needs editing when you change the date. Only add it if the elements come from another time, for example from orbit determination.
- **`orbit.elements`:** Keplerian elements in GCRF. These are currently the MATLAB disposal script's initial orbit: a = 42,464 km, e = 0.005, i = 1°, RAAN and argument of perigee 0. For the start position, give exactly one of:
  - `true_anomaly_deg`, or
  - `start_longitude_deg`: the Earth-fixed longitude at the epoch. The spacecraft stays over the same spot when you change the epoch.
- **`physical_properties`:** mass, plus area and coefficient for drag and for solar radiation pressure (SRP). Both use cannonball (sphere) models.
- **`propagator`:**
  - `force_models`: any of `earth_gravity`, `third_body_sun`, `third_body_moon`, `srp`, `drag`.
  - `gravity_degree` and `gravity_order`: the size of the `earth_gravity` field.
- **`attitude`:**
  - `inertia_kg_m2`: body frame.
  - `initial_quaternion`: `[q0, q1, q2, q3]`, scalar first, rotation GCRF → body.
  - `initial_angular_velocity_deg_s`: body frame.
  - `target_rso`: the name of the RSO that body +Z points at. It must be in the scenario's `rso_file`.
  - `controller`: `natural_frequency_rad_s`, `damping_ratio`.
  - `reaction_wheels`: `max_torque_Nm`, `max_momentum_Nms`, per wheel, with one wheel per body axis. Momentum dumping isn't modeled, so watch `wheel_h*` in `attitude.csv`.
  - `knowledge_error`: `attitude_arcsec` and `rate_deg_s` are 1-sigma errors per axis, redrawn every time step. `seed` makes runs repeatable.
  - `disturbance_torques`: `gravity_gradient` and `srp` (true/false), and `center_of_pressure_offset_m` (body frame, measured from the center of mass).

### RSO.json

The objects propagated with SGP4/SDP4. Each entry has:
- `name`: unique, and used in the outputs.
- `norad_id`, `international_designator` and `object_type` (optional).
- `tle_line1` and `tle_line2`: the exact 69-character lines.

The current three objects are the test set for the orbit starting at 114.7°W, picked from `GEOTLE.json` (TLEs dated Oct 12–13, 2026):

| Object | Type |
|---|---|
| VIASAT_1 | active satellite |
| GOES_3_AKM | rocket body |
| TITAN_3C_TRANSTAGE_DEB_38704 | debris, the pointing target |

### GEOTLE.json

The full GEO-belt catalog, generated by `SpaceTrack/convert_space_track.py` ([steps below](#updating-the-geo-catalog-spacetrack)). Don't edit it by hand: it's overwritten on every refresh.
- **What it covers:** perigee and apogee within ±2,000 km of GEO altitude, inclination under 20°, TLE no older than 30 days at download.
- **Layout:** the same as `RSO.json`. To use it, set `"rso_file": "GEOTLE.json"` in `Scenario.json` (and raise `rso_step_s` first). `RSO.json` stays intact for the test set.
- **Git:** kept out of git by `.gitignore`, because of Space-Track's redistribution terms.

## Outputs

Created by `run.py` in `Outputs/<scenario name>/`. Positions are in GCRF [m, m/s].

| File | Contents |
|---|---|
| `satellite.csv` | Spacecraft position and velocity, every time step |
| `attitude.csv` | Attitude, control, pointing and lighting, every time step (columns below) |
| `rsos.csv` | Every RSO's position and velocity, every `rso_step_s` |
| `rso_summary.csv` | Per RSO: closest approach, % of time visible, % visible and sunlit (closest first) |
| `initial_state.json` | The exact initial state used (the STK verification starts from it) |

`attitude.csv` columns:

| Columns | Meaning |
|---|---|
| `q0`–`q3`, `wx/wy/wz_rad_s` | True attitude (quaternion GCRF → body, scalar first) and body rate |
| `q0_cmd`–`q3_cmd` | Commanded attitude |
| `torque_x/y/z_Nm`, `wheel_hx/hy/hz_Nms` | Wheel torque on the body, and wheel momentum |
| `disturbance_torque_Nm` | Size of the gravity-gradient + solar pressure torque |
| `pointing_error_deg` | Angle between the camera (+Z) and the line of sight to the target |
| `roll_reference` | 0 = roll from the orbit normal, 1 = from the zenith (target near the orbit normal) |
| `target_range_km`, `target_visible`, `target_sunlit` | Distance to the target; 1 if Earth doesn't block the view; 1 if the target is out of Earth's shadow |
| `phase_angle_deg` | Sun–target–spacecraft angle (0 = target fully lit as seen from the spacecraft) |
| `sun_boresight_angle_deg` | Angle between the line of sight and the Sun (small = the camera looks toward the Sun) |
| `spacecraft_lighting` | 0 in Earth's shadow, 1 in full sun |

## Updating the GEO catalog (SpaceTrack/)

Builds `Inputs/GEOTLE.json`, the catalog of everything near the GEO belt. TLEs go stale within days, so refresh before analysis runs:

1. Log in to space-track.org and open this query, which returns the latest TLE for every object whose perigee and apogee are within ±2,000 km of GEO altitude, with inclination under 20°, not decayed, and updated in the last 30 days:
   `https://www.space-track.org/basicspacedata/query/class/gp/PERIAPSIS/33786--37786/APOAPSIS/33786--37786/INCLINATION/<20/DECAY_DATE/null-val/EPOCH/>now-30/orderby/NORAD_CAT_ID/format/json`
2. Save the result as `SpaceTrack/space_track_geo.json`.
3. Run `python SpaceTrack/convert_space_track.py`. It checks every TLE, removes duplicates, gives each object a unique name, and writes `Inputs/GEOTLE.json`.

Neither data file goes into git, because Space-Track's user agreement restricts redistributing its data.

## STK verification

`Verification/stk_verification.py` checks that the simulator is correct by comparing it with STK, an independent tool, using the same inputs. In one run it builds the matching scenario in STK, exports STK's results to `Verification/STK_Results/`, and prints PASS/FAIL.

### How to run it

1. **Run the simulation first:** `python run.py`. The script needs its outputs in `Outputs/<scenario name>/`.
2. **On a Windows machine with STK 12 or 13,** run `python Verification/stk_verification.py`. It attaches to an open STK Desktop or starts one.
   - It needs the STK Python API: `ansys-stk-core` for STK 12.10+/13, or the `agi.stk12` wheel from the STK 12 install.
   - If HPOP can't find a gravity file, pass it explicitly: `--gravity-file "C:\Program Files\AGI\STK 12\STKData\CentralBodies\Earth\EGM2008.grv"`.
   - By default it builds the pointing target plus the 4 RSOs that came closest in the sim. To choose others, use `--rsos NAME NAME ...`.
3. **Read the PASS/FAIL summary** at the end.

To re-run only the comparison on results STK already exported, without opening STK, add `--compare-only`. This works on any machine.

### What STK is set up to match

- **Spacecraft:** HPOP from the exact initial state the sim used, Cartesian in ICRF (= Orekit GCRF), with the force models listed in `Satellite.json`:
  - Earth gravity, same degree and order;
  - Sun and Moon gravity;
  - cannonball solar radiation pressure, same Cr and area/mass, Earth shadow only;
  - drag. It uses a different density model than the sim, but density is essentially zero at GEO.
- **RSOs:** SGP4 from the same TLEs. Only the RSOs STK exported are compared.
- **Attitude:** aligned/constrained, body +Z at the target and body +Y toward the orbit normal. This is the sim's *commanded* attitude. STK has no model of our controller, wheels or sensor errors, so the sim's *true* attitude is checked by confirming its +Z points at the target, using STK's geometry, once the initial slew has settled.

### Checks and limits

The limits are set in the Setup section at the top of `stk_verification.py`.

| Check | Limit | Expect |
|---|---|---|
| Spacecraft position | 100 m | Meters to tens of meters. Kilometers means a force model is on in one tool and off in the other. |
| RSO positions | 1 km | Tens of meters at most. Different SGP4 implementations and TEME conversions. |
| STK +Z vs STK line of sight | 0.01 deg | ~0. Confirms STK's quaternions are being read in our convention. |
| Commanded +Z and +Y vs STK | 0.01 deg + LoS offset | Within the LoS offset |
| True +Z vs STK line of sight, after 600 s | 0.01 deg + LoS offset + 6σ sensor error | ~0.003–0.01 deg (the simulated 10-arcsec attitude knowledge error) |

- **"LoS offset"** is the angle between the sim's and STK's lines of sight to the target. It comes only from their small position differences, and the script prints it. During a close pass, a few meters of position difference becomes a noticeable angle (10 m at 40 km is 0.015 deg).
- **When the roll isn't compared:** while the target is within ~10 deg of the orbit-normal direction, the sim switches its roll reference to the zenith because roll is undefined there. STK keeps using the orbit normal, so the roll check skips those epochs.

### If something fails

- **It has never been run in real STK yet.** It was tested against a scripted fake STK. If the first real run fails, send back the printed output and `Verification/STK_Results/stk_commands.log`.
- **A command was rejected.** `stk_commands.log` lists every STK command and whether STK accepted it. Vector names (`ToTarget`, `Orbit_Normal`, ...) and some options vary between STK versions, so the script tries several variants. If all of them fail, look up the right name in STK and add it to the lists in the script's "Building in STK" section.
- **The convention check fails but the conjugate passes.** The script will say so. It means STK's quaternion convention differs from what's assumed here, not that the sim is wrong.
- **Old results.** Each STK run deletes the old `.e`/`.a` files in `STK_Results/` first, so results from different runs never mix.

## OrbitCalculations

`geo_disposal_hohmann.m` (MATLAB, from the SRR) plans the end-of-mission disposal. It's a two-burn transfer from the operational orbit to a graveyard orbit about 400 km above GEO: about 7.7 m/s total over 12 hours. It also draws a 3D presentation figure. `Inputs/Satellite.json` uses the same starting orbit.

## Status

- [x] Attitude propagation (`Propagation/propagator.py`): pointing, sensors, controller, reaction wheels and disturbance torques
- [x] STK verification script (`Verification/stk_verification.py`). *Next:* run it in STK
- [ ] Replace placeholder spacecraft values with design numbers
- [ ] Run the full `GEOTLE.json` catalog (`rso_file`) after the STK verification
- [ ] Sensor/camera modeling (GitHub issue #2)
- [ ] Reaction wheel momentum dumping
