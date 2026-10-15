import json
import numpy as np

from org.orekit.time import TimeScalesFactory, AbsoluteDate
from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions, Constants
from org.orekit.orbits import KeplerianOrbit, PositionAngleType
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator

from objects import Scenario, RSO, Satellite
from propagator import propagate

if __name__ == '__main__':
    ## Setup frames and time
    EME2000 = FramesFactory.getEME2000()  # J2000 frame
    gcrf = FramesFactory.getGCRF()
    earth_frame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    utc = TimeScalesFactory.getUTC()

    frames = {EME2000: 'EME2000',
              gcrf: 'gcrf',
              earth_frame: 'earth_frame',
              utc: 'utc'}

    # Parse Input
    with open('Inputs/Scenario.json', 'r', encoding='utf8') as file:
        scenario_data = json.load(file)

    scenario_epoch = AbsoluteDate(scenario_data['epoch'], utc)
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
    satellite_epoch = AbsoluteDate(satellite_data['epoch'], utc)
    initial_orbit = KeplerianOrbit(
        semi_major_axis_km * 1000.0,  # semi-major axis in m
        eccentricity,  # eccentricity
        np.deg2rad(inclination_deg),  # inclination
        np.deg2rad(arg_of_perigee_deg),  # argument of perigee
        np.deg2rad(raan_degrees),  # raan
        np.deg2rad(true_anomaly_deg),  # true anomaly
        PositionAngleType.TRUE,  # tells which anomaly type ^
        FramesFactory.getGCRF(),  # frame
        satellite_epoch,  # epoch
        Constants.WGS84_EARTH_MU  # mu
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

    propagate(
        frames=frames,
        sat=sat,
        scenario=scenario,
        rso_propagators=rso_propagators,
        rso_list=rso_list
    )

    print("Sim finished!")