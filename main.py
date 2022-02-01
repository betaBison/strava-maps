"""Create maps of bulk Strava data.

"""

__authors__ = "D. Knowles"
__date__ = "21 Dec 2021"

import os
import math
import gzip
import shutil
import configparser

import fitparse
import matplotlib.pyplot as plt
from gpx_converter import Converter

from lib.plotter import Plotter
from lib.tcx_convert import write_tcx_to_csv
from lib.fit_convert import write_fit_to_csv

def main():

    # read configuration file
    config = configparser.ConfigParser()
    file_dir = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(file_dir,"config","settings.ini")
    config.read(config_path)
    data_path = config["strava"]["data_path"]
    activities_path = os.path.join(data_path,"activities")

    # extract files if necessary
    print("Extracting gz activity files...")
    gz_extract(activities_path)

    # convert gpx files to csv
    print("Converting gpx files to csv...")
    gpx_to_csv(activities_path)

    # convert fit files to csv
    print("Converting fit files to csv...")
    fit_to_csv(activities_path)

    # convert tcx files to csv
    print("Converting tcx files to csv...")
    tcx_to_csv(activities_path)

    # plot maps
    print("Plotting maps...")
    plot_maps(activities_path)


def gz_extract(directory):
    """Extract all files in a directory with gzip.

    Modified from https://gist.github.com/kstreepy/a9800804c21367d5a8bde692318a18f5

    Parameters
    ----------
    directory : string
        Directory path in which file should be extracted.

    """
    full_list = sorted([item for item in os.listdir(directory) if item.endswith(".gz")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, item in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")

        gz_name = os.path.join(directory,item) # get full path of files
        filename = (os.path.basename(gz_name)).rsplit('.',1)[0] #get file name for file withi
        out_path = os.path.join(directory,filename)
        with gzip.open(gz_name,"rb") as f_in, open(out_path,"wb") as f_out:
          shutil.copyfileobj(f_in, f_out)
        os.remove(gz_name) # delete zipped file


def gpx_to_csv(dir):
    """Convert gpx files to csv files.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.

    """
    full_list = sorted([item for item in os.listdir(dir) if item.endswith(".gpx")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, file in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")

        filepath = os.path.join(dir, file)
        outpath = os.path.join(dir, file.rsplit('.',1)[0] + ".csv")
        Converter(input_file=filepath).gpx_to_csv(output_file=outpath)
        os.remove(filepath)

def fit_to_csv(dir):
    """Convert fit files to csv files.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.

    """
    full_list = sorted([item for item in os.listdir(dir) if item.endswith(".fit")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, file in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")
        filepath = os.path.join(dir, file)
        outpath = os.path.join(dir, file.rsplit('.',1)[0] + ".csv")
        fitfile = fitparse.FitFile(filepath,
            data_processor=fitparse.StandardUnitsDataProcessor())
        write_fit_to_csv(fitfile,output_file=outpath)
        os.remove(filepath)

def tcx_to_csv(dir):
    """Convert tcx files to csv files.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.

    """
    full_list = sorted([item for item in os.listdir(dir) if item.endswith(".tcx")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, file in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")
        filepath = os.path.join(dir, file)
        outpath = os.path.join(dir, file.rsplit('.',1)[0] + ".csv")
        write_tcx_to_csv(filepath, outpath)
        os.remove(filepath)

def plot_maps(dir):
    """Map activity data into maps.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.

    """
    p = Plotter(dir)
    # p.stanford()
    # p.bay_area()
    # p.bay_area_laser_elevation_truth()
    p.bay_area_laser_contour_truth()
    plt.show()

if __name__ == "__main__":
    main()
