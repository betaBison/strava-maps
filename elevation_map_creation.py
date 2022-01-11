"""Create merged numpy array of USGS TIFF elevation maps.

"""

__authors__ = "D. Knowles"
__date__ = "27 Dec 2021"

import os
import configparser

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

# sizes in pixels of (width, height)
SIZES = {"small" : (72, 108),
         "medium" : (720, 1080),
         "large" : (7200, 10800)}

def main(img_data, size = "small"):
    """Combine topography maps and save as file.

    Parameters
    ----------
    img_data : dict
        Dictionary that stores topography file names and corresponding
        latitude and longitude data.
    size : string
        Resized shape of each individual elevation map. Options are
        small (default), medium, and large.

    """

    # read configuration file
    config = configparser.ConfigParser()
    file_dir = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(file_dir,"config","settings.ini")
    config.read(config_path)
    data_path = config["usgs"]["data_path"]

    # top row
    im00 = Image.open(os.path.join(data_path, img_data[(0,0)]["filename"]))
    im00 = img_resize(im00, size)
    ar00 = array_convert(im00)
    im01 = Image.open(os.path.join(data_path, img_data[(0,1)]["filename"]))
    im01 = img_resize(im01, size)
    ar01 = array_convert(im01)
    im02 = Image.open(os.path.join(data_path, img_data[(0,2)]["filename"]))
    im02 = img_resize(im02, size)
    ar02 = array_convert(im02)
    ar0 = np.hstack((ar00,ar01,ar02))

    # middle row
    im10 = Image.open(os.path.join(data_path, img_data[(1,0)]["filename"]))
    im10 = img_resize(im10, size)
    ar10 = array_convert(im10)
    im11 = Image.open(os.path.join(data_path, img_data[(1,1)]["filename"]))
    im11 = img_resize(im11, size)
    ar11 = array_convert(im11)
    im12 = Image.open(os.path.join(data_path, img_data[(1,2)]["filename"]))
    im12 = img_resize(im12, size)
    ar12 = array_convert(im12)
    ar1 = np.hstack((ar10,ar11,ar12))

    # bottom row
    ar20 = np.zeros((SIZES[size][1], SIZES[size][0]))
    im21 = Image.open(os.path.join(data_path, img_data[(2,1)]["filename"]))
    im21 = img_resize(im21, size)
    ar21 = array_convert(im21)
    im22 = Image.open(os.path.join(data_path, img_data[(2,2)]["filename"]))
    im22 = img_resize(im22, size)
    ar22 = array_convert(im22)
    ar2 = np.hstack((ar20,ar21,ar22))

    # combine all three rows
    ar_full = np.vstack((ar0,ar1,ar2))

    # save the topography
    output_path = os.path.join(file_dir, "data",
                            "topography_" + size + ".npy")
    np.save(output_path, ar_full)

    plt.figure()
    plt.imshow(ar_full)
    plt.show()


def img_resize(im, size):
    """Resize images appropriately.

    Parameters
    ----------
    im : pillow image
        Image to be resized
    size : string
        Resized shape of each individual elevation map. Options are
        small (default), medium, and large.

    Returns
    -------
    im_resized : pillow image
        resized image


    """
    im_resized = im.resize(SIZES[size])

    return im_resized

def array_convert(im):
    """Convert pillow image to numpy array.

    Parameters
    ----------
    im : pillow image
        Image to be resized

    Returns
    -------
    img_array : np.ndarray
        numpy array

    """
    img_array = np.array(im)
    img_array = np.clip(img_array, 0.0, np.inf)

    return img_array


if __name__ == "__main__":
    img_data = {(1,1):{"filename":"USGS_13_n38w123_20210617.tif","minX":-123.00055555579371,"maxX":-121.99944444360563,"minY":36.999444443707034,"maxY":38.00055555589506},
                (1,2):{"filename":"USGS_13_n38w122_20210617.tif","minX":-122.00055555599357,"maxX":-120.99944444380549,"minY":36.999444443707034,"maxY":38.00055555589506},
                (2,2):{"filename":"USGS_13_n37w122_20210617.tif","minX":-122.00055555599357,"maxX":-120.99944444380549,"minY":35.999444443906840,"maxY":37.00055555609492},
                (2,1):{"filename":"USGS_13_n37w123_20210617.tif","minX":-123.00055555579371,"maxX":-121.99944444360563,"minY":35.999444443906840,"maxY":37.00055555609492},
                (0,0):{"filename":"USGS_13_n39w124_20200115.tif","minX":-124.00055555649300,"maxX":-122.99944444340600,"minY":37.999444443507200,"maxY":39.00055555659460},
                (0,1):{"filename":"USGS_13_n39w123_20200317.tif","minX":-123.000555555794,"maxX":-121.999444443606,"minY":37.9994444435072,"maxY":39.0005555565946},
                (0,2):{"filename":"USGS_13_n39w122_20200317.tif","minX":-122.000555555994,"maxX":-120.999444443805,"minY":37.9994444435072,"maxY":39.0005555565946},
                (1,0):{"filename":"USGS_13_n38w124_20200115.tif","minX":-124.000555556493,"maxX":-122.999444443406,"minY":36.999444443707,"maxY":38.0005555558951},
                }

    main(img_data, "medium")
