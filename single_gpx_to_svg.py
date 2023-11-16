"""Convert a single GPX file to svg.


"""


__authors__ = "D. Knowles"
__date__ = "15 Nov 2023"

import pandas as pd
import matplotlib.pyplot as plt
from gpx_converter import Converter

def main():

    filepath = "/home/<path>/<filename>.gpx"

    # convert gps to csv
    outpath = filepath[:-3] + "csv"
    Converter(input_file=filepath).gpx_to_csv(output_file=outpath)

    print("inputing csv data")
    df = pd.read_csv(outpath)

    fig,ax = plt.subplots(1,1)

    plt.plot(df["longitude"],df["latitude"])

    ax.axis('equal')

    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.set_axis_off()
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
        hspace = 0, wspace = 0)
    plt.margins(0,0)
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    fig.savefig(filepath[:-3]+"svg",
                format="svg",
                dpi=300.,
                pad_inches=0,
                bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
