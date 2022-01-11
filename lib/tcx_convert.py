"""Convert tcx files to csv files.

Copied from https://github.com/coreysiegel/tcx-gpx-csv/blob/master/tcx2csv.py

"""

__authors__ = "C. Siegel, D. Knowles"
__date__ = "23 Dec 2021"

import re
import sys
import pytz
import string
from datetime import datetime
import xml.etree.ElementTree as et

UTC = pytz.UTC

# takes in a TCX file and outputs a CSV file
def write_tcx_to_csv(input, output):
    """Converts .tcx files to .csv files.

    Parameters
    ----------
    input : string
        Input filepath.
    output : string
        Output filepath.

    """
    
    strip_whitespace(input)

    tree = et.parse(input)
    root = tree.getroot()
    m=re.match(r'^({.*})', root.tag)
    if m:
        ns=m.group(1)
    else:
        ns=''
    if root.tag!=ns+'TrainingCenterDatabase':
        print('Unknown root found: '+root.tag)
        return
    activities=root.find(ns+'Activities')
    if not activities:
        print('Unable to find Activities under root')
        return
    activity=activities.find(ns+'Activity')
    if not activity:
        print('Unable to find Activity under Activities')
        return
    columnsEstablished=False
    for lap in activity.iter(ns+'Lap'):
        if columnsEstablished:
            fout.write('New Lap\n')
        for track in lap.iter(ns+'Track'):
            if columnsEstablished:
                fout.write('New Track\n')
            for trackpoint in track.iter(ns+'Trackpoint'):
                try:
                    string_time=trackpoint.find(ns+'Time').text.strip()
                    dtime = datetime.strptime(string_time, "%Y-%m-%dT%H:%M:%S.%f+00:00")
                    dtime = UTC.localize(dtime).astimezone(UTC)
                    time = str(dtime.replace(microsecond=0).isoformat(' '))
                except:
                    time = ''
                try:
                    latitude=trackpoint.find(ns+'Position').find(ns+'LatitudeDegrees').text.strip()
                    latitude=str(round(float(latitude),6))
                except:
                    latitude=''
                try:
                    longitude=trackpoint.find(ns+'Position').find(ns+'LongitudeDegrees').text.strip()
                    longitude=str(round(float(longitude),6))
                except:
                    longitude=''
                try:
                    altitude=trackpoint.find(ns+'AltitudeMeters').text.strip()
                    altitude=str(round(float(altitude),1))
                except:
                    altitude=''
                if not columnsEstablished:
                    fout=open(output, 'w')
                    fout.write(','.join(('time','latitude','longitude','altitude'))+'\n')
                    columnsEstablished=True
                fout.write(','.join((time,latitude,longitude,altitude))+'\n')

    fout.close()


def strip_whitespace(filename):
    """Strips whitespace from the file.

    For some reason, my .tcx files had some extra whitespace in the
    beginning of the files that caused errors in the tcx conversion, so
    before converting, this removed the whitespace at the beginning of
    the file.

    Parameters
    ----------
    filename : string
        Name of the file that should be stripped of whitespace.

    """
    file = open(filename, "r")
    contents = file.read()
    new_contents = contents.strip()
    file.close()

    file = open(filename, "w")
    file.write(new_contents)
    file.close()
