"""Plot strava map data

"""

__authors__ = "D. Knowles"
__date__ = "21 Dec 2021"

import os
import sys
import math
import time
# append <path>/strave-maps/ to path
sys.path.append(os.path.dirname(
                os.path.realpath(__file__)))

import numpy as np
import pandas as pd
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import LineCollection

# import specific parameters
# from config.contours import Tristan as Params
from config.contours import Derek as Params

class Plotter():
    def __init__(self, dir):
        """Plot fun maps from Strava data.

        Parameters
        ----------
        dir : string
            Directory path in which file should be extracted.

        """

        self.P = Params()

        file_dir = os.path.dirname(os.path.realpath(__file__))
        self.repo_dir = os.path.join(file_dir, "..")

        self.df = []

        self.plot_time = str(time.time())[:10]

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

        fig = plt.figure(figsize=self.P.figsize)

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
        ctopo = ndimage.gaussian_filter(ctopo, self.P.blur)
        plt.contourf(x,y,ctopo, self.P.levels,
                     colors=self.P.colors)

        for df in self.df:
            plt.plot(df["longitude"].to_numpy(), df["latitude"].to_numpy(),"k")

        plt.xlim(self.P.xlim)
        plt.ylim(self.P.ylim)

        # remove axis ticks and labels
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

        fig.savefig(os.path.join(self.repo_dir,"plots",
                    self.P.name + "-proof-" + self.plot_time \
                    + ".png"),
                    format="png",
                    dpi=300.,
                    bbox_inches="tight")
        plt.close(fig)

    def bay_area_laser_elevation_truth(self):
        """Laser contours for the entire Bay area."""

        levels = self.P.levels

        fig = plt.figure(figsize=self.P.figsize)
        ax = plt.subplot(1,1,1)
        # plt.figure(figsize=(8.,10.))

        topo = np.load(os.path.join(self.repo_dir,"data",
                                    "topography_medium.npy"))
        left = -124.000555556493
        right = -120.99944444380549
        bottom = 35.999444443906840
        top = 39.0005555565946

        y_range = np.linspace(bottom,top,topo.shape[0])
        x_range = np.linspace(left,right,topo.shape[1])
        x,y = np.meshgrid(x_range,y_range)
        ctopo = np.flip(topo,axis=0)
        # ctopo = ndimage.gaussian_filter(ctopo, 3)


        for dd, df in enumerate(self.df):

            if len(df) == 0:
                continue

            # update colors
            df["color"] = "k"
            for ll in range(len(levels) - 1):
                df.loc[(df["altitude"] > levels[ll]) & (df["altitude"] <= levels[ll + 1]), "color"] = "C" + str(ll)
                # print(len(sub_df),levels[ll],levels[ll+1])

            latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
                                  df["latitude"].to_numpy().reshape(-1,1)))

            sequence = np.zeros((latlon.shape[0] - 1, 2, 2))
            sequence[:,0,:] = latlon[:-1,:]
            sequence[:,1,:] = latlon[1:,:]
            print(sequence.shape)
            lc = LineCollection(sequence,
                                colors = df["color"].to_numpy())

            lc.set_linewidth(2)
            line = ax.add_collection(lc)

        plt.contour(x,y,ctopo, levels = levels,
                    colors="k")

        plt.xlim(self.P.xlim)
        plt.ylim(self.P.ylim)

        # remove axis ticks and labels
        ax = plt.gca()
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)

        fig.savefig(os.path.join(self.repo_dir,"plots",
                    self.P.name + "-elevation-" + self.plot_time\
                     + ".png"),
                    format="png",
                    dpi=300.,
                    bbox_inches="tight")
        plt.close(fig)

    def bay_area_laser_contours(self):
        """Laser contours for the entire Bay area.

        """
        levels = self.P.levels
        level_indexes = [3, 2, 1]
        # self.df = self.df[1000:1400]

        for dd, df in enumerate(self.df):
            df["color"] = "orange"

        proof_fig = plt.figure(figsize=self.P.figsize)
        proof_ax = plt.subplot(1,1,1)

        contoured_paths = []

        for level_index in level_indexes:
        # for level_index in range(1):

            fig = plt.figure(figsize=self.P.figsize)
            ax = plt.subplot(1,1,1)
            # plt.figure(figsize=(8.,10.))

            topo = np.load(os.path.join(self.repo_dir,"data",
                                        "topography_medium.npy"))
            left = -124.000555556493
            right = -120.99944444380549
            bottom = 35.999444443906840
            top = 39.0005555565946

            y_range = np.linspace(bottom,top,topo.shape[0])
            x_range = np.linspace(left,right,topo.shape[1])
            x,y = np.meshgrid(x_range,y_range)
            ctopo = np.flip(topo,axis=0)
            ctopo = ndimage.gaussian_filter(ctopo, self.P.blur)

            contour_colors = [[0.,0.,0.,0.]]*len(levels)
            # contour_colors[level_index + 1] = "k"

            CS = plt.contour(x,y,ctopo, levels = levels,
                        colors=contour_colors)

            # remove small collections unless they're in the white list
            filtered_paths = []
            collection = CS.collections[level_index]
            for path in collection.get_paths():
                if len(path.vertices) > 500 and level_index == 1:
                    # add extra vertices for plotting for largest path
                    new_vertices = path.vertices.copy()
                    top_left = new_vertices[0:1,:]
                    top_left[0,1] += 1
                    bottom_right = new_vertices[-1,:].reshape(1, -1)
                    bottom_right[0,0] += 1
                    top_right = np.array([[bottom_right[0,0], top_left[0,1]]])
                    new_vertices = np.insert(new_vertices, [0],
                                             top_left, axis = 0)
                    new_vertices = np.append(new_vertices,
                                             bottom_right, axis = 0)
                    new_vertices = np.append(new_vertices,
                                             top_right, axis = 0)
                    new_path = Path(new_vertices)
                    path = new_path
                if len(path.vertices) < self.P.contour_min:
                    ignore = True
                    if level_index == 1:
                        center_location = np.mean(path.vertices,axis=0)[::-1]
                        white_list = self.P.wlc
                        for ww in range(white_list.shape[0]):
                            diff = self.haversine(center_location, white_list[ww])
                            if diff < 5000:
                                ignore = False
                                break
                    if ignore:
                        continue


                filtered_paths.append(path)
                patch = patches.PathPatch(path, edgecolor='k', facecolor='w', lw=2)
                ax.add_patch(patch)

            contoured_paths.insert(0,filtered_paths)

            # percent10 = int(len(self.df)/10.)
            # for dd, df in enumerate(self.df):
            #     # print("dd:",dd,"/",len(self.df))
            #     if (dd+1) % percent10 == 0:
            #         print(math.ceil(100*float(dd)/len(self.df)),"% of activities plotted")
            #
            #     if len(df) == 0:
            #         continue
            #
            #     # # update colors
            #     # df["color"] = "orange"
            #     # for ll in range(len(levels) - 1):
            #     #     df.loc[(df["altitude"] > levels[ll]) & (df["altitude"] <= levels[ll + 1]), "color"] = "C" + str(ll)
            #     #     # print(len(sub_df),levels[ll],levels[ll+1])
            #
            #     latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
            #                           df["latitude"].to_numpy().reshape(-1,1)))
            #
            #     sequence = np.zeros((latlon.shape[0] - 1, 2, 2))
            #     sequence[:,0,:] = latlon[:-1,:]
            #     sequence[:,1,:] = latlon[1:,:]
            #
            #     # for ss in range(sequence.shape[0]):
            #     #     for cc, collection in enumerate(CS.collections):
            #     #         for path in collection.get_paths():
            #     #             if np.any(path.contains_points(sequence[ss,:,:])):
            #     #                 df.loc[ss,"color"] = "C" + str(cc)
            #
            #     # color_options = ["m","b","g","r","m"]
            #
            #     for path in filtered_paths:
            #         mask = path.contains_points(latlon)
            #         # if np.any(mask):
            #         #     print(mask.shape)
            #         #     print(mask)
            #         df.loc[mask,"color"] = "b"
            #
            #     lc = LineCollection(sequence,
            #                         colors = df["color"].to_numpy())
            #
            #     lc.set_linewidth(1)
            #     line = ax.add_collection(lc)


            plt.xlim(self.P.xlim)
            plt.ylim(self.P.ylim)

            # remove axis ticks and labels
            ax = plt.gca()
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)

            fig.savefig(os.path.join(self.repo_dir,"plots",
                        self.P.name + "-contour-" \
                      + self.plot_time \
                      + "-" + str(level_index) + ".png"),
                        format="png",
                        dpi=300.,
                        bbox_inches="tight")
            plt.close(fig)

        print("Done with all contour vectors.\n")


        ################################################################
        # visual contour proof
        ################################################################
        background_vertices = np.array([[self.P.xlim[0],self.P.ylim[0]],
                                        [self.P.xlim[0],self.P.ylim[1]],
                                        [self.P.xlim[1],self.P.ylim[1]],
                                        [self.P.xlim[1],self.P.ylim[0]],
                                        ])
        background_path = Path(background_vertices)
        background = patches.PathPatch(background_path, edgecolor='k',
                                        facecolor=self.P.colors[0], lw=2)
        proof_ax.add_patch(background)

        for level_index in range(3):
            """Add new paths."""

            hole_patches = []
            for pp, path in enumerate(contoured_paths[level_index]):

                for cp, check_path in enumerate(contoured_paths[level_index]):
                    contained_measure = sum(check_path.contains_points(path.vertices))/float(path.vertices.shape[0])
                    # print(check_path.contains_points(path.vertices))
                    # print(contained_measure)
                    if contained_measure > self.P.contained_measure and pp != cp:
                        color = self.P.colors[level_index]
                        #
                        # print(pp,cp,"hey",color,level_index)
                        # ax_new = plt.gca()
                        # # proof fig
                        # ax_new.set_xlim(self.P.xlim)
                        # ax_new.set_ylim(self.P.ylim)
                        #
                        # # remove axis ticks and labels
                        # ax_new.axes.xaxis.set_visible(False)
                        # ax_new.axes.yaxis.set_visible(False)
                        new_patch = patches.PathPatch(path, edgecolor="k",
                                                        facecolor=color, lw=0)
                        hole_patches.append(new_patch)
                        # ax_new.add_patch(new_patch)
                        # plt.show()

                    else:
                        color = self.P.colors[level_index + 1]

                    proof_patch = patches.PathPatch(path, edgecolor="k",
                                                    facecolor = color, lw=0)
                    proof_ax.add_patch(proof_patch)


            for hole_patch in hole_patches:
                proof_ax.add_patch(hole_patch)

            print("added proof contours for level",level_index+1,"/3")


        percent10 = int(len(self.df)/10.)
        for dd, df in enumerate(self.df):
            # print("dd:",dd,"/",len(self.df))
            if (dd+1) % percent10 == 0:
                print(math.ceil(100*float(dd)/len(self.df)),"% of activities plotted")

            if len(df) == 0:
                continue

            if len(df) > self.P.path_min:
                plt.plot(df["longitude"].to_numpy().reshape(-1,1),
                         df["latitude"].to_numpy().reshape(-1,1),
                         "k",linewidth = 0.5)
            # else:
            #     plt.plot(df["longitude"].to_numpy().reshape(-1,1),
            #     df["latitude"].to_numpy().reshape(-1,1),
            #     "b",linewidth = 0.3)


        # proof fig
        proof_ax.set_xlim(self.P.xlim)
        proof_ax.set_ylim(self.P.ylim)

        # remove axis ticks and labels
        proof_ax = plt.gca()
        proof_ax.axes.xaxis.set_visible(False)
        proof_ax.axes.yaxis.set_visible(False)

        proof_fig.savefig(os.path.join(self.repo_dir,"plots",
                    self.P.name + "-contour-proof-" \
                  + self.plot_time + ".png"),
                    format="png",
                    dpi=300.,
                    bbox_inches="tight")
        plt.close(proof_fig)

    def haversine(self, coord1, coord2):
        """
        Taken from: https://janakiev.com/blog/gps-points-distance-python/
        """
        R = 6372800  # Earth radius in meters
        lat1, lon1 = coord1
        lat2, lon2 = coord2

        phi1, phi2 = math.radians(lat1), math.radians(lat2)
        dphi       = math.radians(lat2 - lat1)
        dlambda    = math.radians(lon2 - lon1)

        a = math.sin(dphi/2)**2 + \
            math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2

        return 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))
