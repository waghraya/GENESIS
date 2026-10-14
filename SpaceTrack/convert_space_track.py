# Converts a raw Space-Track "gp" JSON download into Inputs/GEOTLE.json (same layout as Inputs/RSO.json, plus object
# type, radar size class and element set epoch). Every TLE is checked; bad entries and duplicates are dropped.
#     python SpaceTrack/convert_space_track.py [download.json] [output.json]
# Defaults: SpaceTrack/space_track_geo.json -> Inputs/GEOTLE.json. The README has the Space-Track query.
import collections
import json
import os
import re
import sys

SPACE_TRACK_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(os.path.dirname(SPACE_TRACK_DIR), 'Inputs')


# TLE checksum: mod-10 sum of the digits in columns 1-68, with '-' counting as 1 and everything else 0.
# Column 69 holds the expected value, so a mistyped or truncated line almost always fails it.
def tle_checksum(line):
    return sum(int(c) if c.isdigit() else (1 if c == '-' else 0) for c in line[:68]) % 10


# Reason the TLE pair is malformed, or None if it is fine
def tle_problem(line1, line2):
    if len(line1) != 69 or len(line2) != 69:
        return 'line length is not 69'
    if not (line1.startswith('1 ') and line2.startswith('2 ')):
        return 'line numbers'
    if not (line1[68].isdigit() and line2[68].isdigit()) \
            or int(line1[68]) != tle_checksum(line1) or int(line2[68]) != tle_checksum(line2):
        return 'checksum'
    if line1[2:7] != line2[2:7]:
        return 'catalog numbers differ between lines'
    return None


# Names label the simulation outputs and become STK object names: keep letters, digits, '-' and '_' only
def clean_name(name):
    return re.sub(r'_+', '_', re.sub(r'[^A-Za-z0-9_-]', '_', name)).strip('_')


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SPACE_TRACK_DIR, 'space_track_geo.json')
    output_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(INPUT_DIR, 'GEOTLE.json')
    with open(input_path, encoding='utf-8') as file:
        records = json.load(file)

    # Keep the newest valid element set per object
    newest, skipped = {}, collections.Counter()
    for record in records:
        problem = tle_problem(record['TLE_LINE1'], record['TLE_LINE2'])
        if problem:
            skipped[problem] += 1
            continue
        norad_id = int(record['NORAD_CAT_ID'])
        if norad_id not in newest or record['EPOCH'] > newest[norad_id]['EPOCH']:
            newest[norad_id] = record
    kept = [newest[norad_id] for norad_id in sorted(newest)]

    # Names must be unique (they identify each RSO in the outputs): number the repeated ones with their NORAD ID
    name_counts = collections.Counter(clean_name(record['OBJECT_NAME']) for record in kept)
    rsos = [{
        'name': clean_name(record['OBJECT_NAME']) if name_counts[clean_name(record['OBJECT_NAME'])] == 1
        else f"{clean_name(record['OBJECT_NAME'])}_{record['NORAD_CAT_ID']}",
        'norad_id': record['NORAD_CAT_ID'],
        'international_designator': record['OBJECT_ID'],
        'object_type': record['OBJECT_TYPE'],
        'rcs_size': record['RCS_SIZE'],
        'epoch': record['EPOCH'],
        'tle_line1': record['TLE_LINE1'],
        'tle_line2': record['TLE_LINE2'],
    } for record in kept]

    with open(output_path, 'w', encoding='utf-8') as file:
        json.dump({'source': 'space-track.org gp class', 'created_utc': max(record['CREATION_DATE'] for record in kept),
                   'count': len(rsos), 'rsos': rsos}, file, indent=2)

    print(f'{len(records)} records read, {len(rsos)} objects written to {output_path}')
    print(f'  skipped: {dict(skipped) if skipped else "none"}, duplicates removed: {len(records) - sum(skipped.values()) - len(rsos)}')
    print(f'  by type: {dict(collections.Counter(rso["object_type"] for rso in rsos))}')
    print(f'  element set epochs: {min(rso["epoch"] for rso in rsos)} to {max(rso["epoch"] for rso in rsos)}')


if __name__ == '__main__':
    main()
