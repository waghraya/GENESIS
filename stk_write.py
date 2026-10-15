from java.io import PrintStream
from org.orekit.time import AbsoluteDate, TimeScalesFactory

def write_stk_e(ephemeris_provider, file_path: str, end_date: AbsoluteDate, start_date: AbsoluteDate, time_step_seconds: float):
    ps = PrintStream(file_path)

    current_date = start_date
    # TODO: make dependent on input file values
    # Ephemeris header
    ps.println("stk.v.12.0")  # Adjust version string if needed
    ps.println("BEGIN Ephemeris")
    ps.println("NumberOfEphemerisPoints " + str(int((end_date.durationFrom(start_date)) / time_step_seconds) + 1))
    ps.println("InterpolationMethod     Lrange")
    ps.println("InterpolationOrder      5")
    ps.println("CentralBody             Earth")
    ps.println("CoordinateSystem        ICRF")  # Adjust based on your Orekit output frame (e.g., GCRF/EME2000)
    ps.println("DistanceUnit            Kilometers")
    ps.println("TimeFormat              UTCG")  # Using Gregorian UTC (e.g., 22 May 2026 12:00:00.000)
    ps.println("EphemerisTimePosVel")
    ps.println()
    ps.println("Data")
    # Iterate and write data points
    while current_date.compareTo(end_date) <= 0:

    # Write footer


    ps.close()