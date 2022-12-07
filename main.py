"""Create maps of bulk Strava data.

"""

__authors__ = "D. Knowles"
__date__ = "21 Dec 2021"

import os
import math
import gzip
import shutil

import fitparse
import matplotlib.pyplot as plt
from gpx_converter import Converter

from lib.plotter import Plotter
from lib.tcx_convert import write_tcx_to_csv
from lib.fit_convert import write_fit_to_csv
# from config.settings import Template as Settings
from config.settings import Derek as Settings
# from config.settings import Tristan as Settings

def main():
    S = Settings()
    data_path = S.data_path
    activities_path = os.path.join(data_path,"activities")
    activities_output = os.path.join(data_path,"activities_csv")
    make_dir(activities_output)
    print("path:",activities_path)

    # extract files if necessary
    print("Extracting gz activity files...")
    gz_extract(activities_path)

    # convert gpx files to csv
    print("Converting gpx files to csv...")
    gpx_to_csv(activities_path, activities_output)

    # convert fit files to csv
    print("Converting fit files to csv...")
    fit_to_csv(activities_path, activities_output)

    # convert tcx files to csv
    print("Converting tcx files to csv...")
    tcx_to_csv(activities_path, activities_output)

    # plot maps
    print("Plotting maps...")
    plot_maps(S)


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

def gpx_to_csv(dir, activities_output):
    """Convert gpx files to csv files.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.
    activities_output : string
        Directory path in which file should be saved.

    """
    full_list = sorted([item for item in os.listdir(dir) if item.endswith(".gpx")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, file in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")

        filepath = os.path.join(dir, file)
        outpath = os.path.join(activities_output, file.rsplit('.',1)[0] + ".csv")
        Converter(input_file=filepath).gpx_to_csv(output_file=outpath)

def fit_to_csv(dir, activities_output):
    """Convert fit files to csv files.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.
    activities_output : string
        Directory path in which file should be saved.

    """
    full_list = sorted([item for item in os.listdir(dir) if item.endswith(".fit")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, file in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")
        filepath = os.path.join(dir, file)
        outpath = os.path.join(activities_output, file.rsplit('.',1)[0] + ".csv")
        fitfile = fitparse.FitFile(filepath,
            data_processor=fitparse.StandardUnitsDataProcessor())
        write_fit_to_csv(fitfile,output_file=outpath)

def tcx_to_csv(dir, activities_output):
    """Convert tcx files to csv files.

    Parameters
    ----------
    dir : string
        Directory path in which file should be extracted.
    activities_output : string
        Directory path in which file should be saved.

    """
    full_list = sorted([item for item in os.listdir(dir) if item.endswith(".tcx")])
    percent10 = math.ceil(len(full_list)/10.)

    for ii, file in enumerate(full_list):
        if ii % percent10 == 0:
            print(int(100*float(ii)/len(full_list)),"% complete")
        filepath = os.path.join(dir, file)
        outpath = os.path.join(activities_output, file.rsplit('.',1)[0] + ".csv")
        write_tcx_to_csv(filepath, outpath)

def plot_maps(S):
    """Map activity data into maps.

    Parameters
    ----------
    S : settings class instance
        User settings class instance.

    """
    p = Plotter(S)
    p.laser_contours()

def make_dir(directory): # pragma: no cover
    """Create a file directory if it doesn't yet exist.
    Parameters
    ----------
    directory : string
        Filepath of directory to create if it does not exist.
    """

    # create directory if it doesn't yet exist
    if not os.path.isdir(directory):
        try:
            os.makedirs(directory)
        except OSError as error:
            raise OSError("Unable to create directory " + directory) from error

if __name__ == "__main__":
    main()
