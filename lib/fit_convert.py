"""Convert .fit files to .csv

Modified from Max Candocia's work at
https://github.com/mcandocia/fit_processing
and
https://github.com/mcandocia/fit_processing/blob/master/convert_fit_to_csv.py

"""

__authors__ = "M. Candocia, D. Knowles"
__date__ = "23 Dec 2021"

import os
import csv
import pytz
import fitparse

allowed_fields = ['timestamp','position_lat','position_long', 'distance',
'enhanced_altitude', 'altitude','enhanced_speed',
                 'speed', 'heart_rate','cadence','fractional_cadence']
required_fields = ['timestamp', 'position_lat', 'position_long']
label_names = ["time", "latitude", "longitude", "altitude"]

UTC = pytz.UTC

def main():
    fitfile = fitparse.FitFile("XXX.fit")

    print('converting')
    write_fit_to_csv(fitfile)
    print('finished conversions')


def write_fit_to_csv(fitfile, output_file='test_output.csv'):
    """Convert .fit files to .csv files.

    Parameters
    ----------
    fitfile : fitfile
        The fitfile to be translated.
    output_file : string
        The name of the output csv file.

    """

    csv_data = []
    for rr, record in enumerate(fitfile.get_messages("record")):
        # Records can contain multiple pieces of data (ex: timestamp, latitude, longitude, etc)
        mdata = {}
        for data in record:
            if data.name in required_fields:
                if data.value is None:
                    continue
                if data.name=='timestamp':
                    timestamp = UTC.localize(data.value).astimezone(UTC)
                    mdata[data.name] = timestamp.strftime("%Y:%m:%d:%H:%M:%S")
                elif data.name in ("position_lat", "position_long"):
                    if data.units == "semicircles":
                        degrees = float(data.value) * (180. / 2**31)
                        mdata[data.name] = round(degrees, 6)
                    elif data.units == "deg":
                        mdata[data.name] = data.value
                    else:
                        print("unknown lat/long units:",data.units)
                else:
                    mdata[data.name] = data.value

        write_data = True
        for key in required_fields:
            if key not in mdata:
                write_data = False
        if write_data:
            csv_data.append(mdata)

    #write to csv
    with open(output_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(label_names)
        for entry in csv_data:
            writer.writerow([ str(entry.get(k, '')) for k in required_fields])

if __name__=='__main__':
    main()
