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
required_fields = ['timestamp', 'position_lat', 'position_long', 'altitude']
label_names = ["time", "latitude", "longitude", "altitude"]

UTC = pytz.UTC

def main():
    fitfile = fitparse.FitFile("XXX.fit",
        data_processor=fitparse.StandardUnitsDataProcessor())

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



    messages = fitfile.messages
    data = []
    for m in messages:
        skip=False
        if not hasattr(m, 'fields'):
            continue
        fields = m.fields
        #check for important data types
        mdata = {}
        for field in fields:
            if field.name in allowed_fields:
                if field.value == None:
                    skip=True
                elif field.name=='timestamp':
                    mdata[field.name] = UTC.localize(field.value).astimezone(UTC)
                elif field.name=="position_lat" or field.name == "position_long":
                    mdata[field.name] = round(field.value, 6)
                else:
                    mdata[field.name] = field.value
        for rf in required_fields:
            if rf not in mdata:
                skip=True
        if not skip:
            data.append(mdata)
    #write to csv
    with open(output_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(label_names)
        for entry in data:
            writer.writerow([ str(entry.get(k, '')) for k in required_fields])
    # print('wrote %s' % output_file)

if __name__=='__main__':
    main()
