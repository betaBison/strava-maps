"""Plot strava map data

"""

__authors__ = "D. Knowles"
__date__ = "21 Dec 2021"

import os
import math

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt

class Plotter():
    def __init__(self, dir):
        """Plot fun maps from Strava data.

        Parameters
        ----------
        dir : string
            Directory path in which file should be extracted.

        """

        file_dir = os.path.dirname(os.path.realpath(__file__))
        self.repo_dir = os.path.join(file_dir, "..")

        self.df = []

        full_list = sorted([item for item in os.listdir(dir) if item.endswith(".csv")])
        percent10 = int(len(full_list)/10.)

        print("inputing csv data")
        for ii, file in enumerate(full_list):
            if (ii+1) % percent10 == 0:
                print(math.ceil(100*float(ii)/len(full_list)),"% complete")
            filepath = os.path.join(dir, file)
            new_df = pd.read_csv(filepath)
            self.df.append(new_df)


    def stanford(self):
        """Plot area near Stanford."""

        for df in self.df:
            plt.plot(df["longitude"].to_numpy(), df["latitude"].to_numpy())

        plt.xlim([-122.2017295427098,-122.13621255910668])
        plt.ylim([37.40984228647541, 37.44176440070531])

        # remove axis ticks and labels
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)


    def bay_area(self):
        """Plot of entire Bay area."""

        # plt.figure(figsize=(4.8,7.2))
        # plt.figure(figsize=(9.6,14.4))
        # plt.figure(figsize=(24.,36.))
        fig = plt.figure(figsize=(16.,20.))
        # plt.figure(figsize=(8.,10.))

        topo = np.load(os.path.join(self.repo_dir,"data",
                                    "topography_medium.npy"))
        left = -124.000555556493
        right = -120.99944444380549
        bottom = 35.999444443906840
        top = 39.0005555565946
        # plt.imshow(topo, extent = (left, right, bottom, top))

        y_range = np.linspace(bottom,top,topo.shape[0])
        x_range = np.linspace(left,right,topo.shape[1])
        x,y = np.meshgrid(x_range,y_range)
        ctopo = np.flip(topo,axis=0)
        ctopo = ndimage.gaussian_filter(ctopo, 3)
        plt.contourf(x,y,ctopo, [-1000, 2., 150., 350, 5000],
                     colors=["#afeeee","#7B3F00","#D27D2D","#C19A6B"])

        for df in self.df:
            plt.plot(df["longitude"].to_numpy(), df["latitude"].to_numpy(),"k")

        plt.xlim([-122.71923426539561,-121.7753358765508])
        plt.ylim([36.94154503938612, 38.04094280129413])

        # remove axis ticks and labels
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

        fig.savefig("bay_area.png",
                    format="png",
                    bbox_inches="tight")
        plt.close(fig)
