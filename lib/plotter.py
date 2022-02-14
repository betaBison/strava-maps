"""Plot strava map data.

"""

__authors__ = "D. Knowles"
__date__ = "21 Dec 2021"

import os
import sys
import math
import time
import copy
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

class Plotter():
    def __init__(self, S):
        """Plot fun maps from Strava data.

        Parameters
        ----------
        S : settings class instance
            User settings class instance.

        """

        self.S = S
        data_path = self.S.data_path
        activities_path = os.path.join(data_path,"activities")

        file_dir = os.path.dirname(os.path.realpath(__file__))
        self.repo_dir = os.path.join(file_dir, "..")

        self.df = []

        self.plot_time = str(time.time())[:10]

        full_list = sorted([item for item in os.listdir(activities_path) if item.endswith(".csv")])
        percent10 = int(len(full_list)/10.)

        print("inputing csv data")
        for ii, file in enumerate(full_list):
            if (ii+1) % percent10 == 0:
                print(math.ceil(100*float(ii)/len(full_list)),"% complete")
            filepath = os.path.join(activities_path, file)
            new_df = pd.read_csv(filepath)
            if len(new_df) > 0:
                self.df.append(new_df)

    def get_contours(self):
        """Get the contours from the topography input map.

        Adds contours to a dummy figure that is immediately closed.

        Returns
        -------
        contours : QuadContourSet
            Contour object see, plt.contour.QuadContourSet for
            more documentation. This is the return of plt.contour()

        """

        # load topography from previously saved numpy array in the
        # data/ directory
        topo = np.load(os.path.join(self.repo_dir,"data",
                                    self.S.topo_name))
        # outside coordinates of the topography map
        x_range = np.linspace(self.S.topo_xlim[0], self.S.topo_xlim[1],
                              topo.shape[1])
        y_range = np.linspace(self.S.topo_ylim[0], self.S.topo_ylim[1],
                              topo.shape[0])
        x,y = np.meshgrid(x_range,y_range)
        # flip to change from array to image
        ctopo = np.flip(topo,axis=0)
        # gaussian blur topography to get smoother effect
        ctopo = ndimage.gaussian_filter(ctopo, self.S.blur)

        # add contours to temporary figure that's closed immediately
        temp_fig = plt.figure(figsize=self.S.figsize)
        contours = plt.contour(x, y, ctopo, levels = self.S.levels, colors="red")
        self.save_and_close(temp_fig, "temp")

        return contours

    def get_contour_paths(self):
        """Get contour paths for all levels.

        Returns
        -------
        contour_paths : list
            List of lists of all contour paths.
            Order of [[1st level paths],[2nd level paths],[etc...]].

        """

        contour_paths = []
        contours = self.get_contours()

        for level_index in range(1,len(self.S.levels)-1):
            # remove small collections unless they're in the white list
            filtered_paths = []
            collection = contours.collections[level_index]
            for path in collection.get_paths():
                if len(path.vertices) > 500 and level_index == 1:
                    # add extra vertices for plotting the largest path
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
                if len(path.vertices) < self.S.contour_min:
                    ignore = True
                    if level_index == 1:
                        center_location = np.mean(path.vertices,axis=0)[::-1]
                        white_list = self.S.wlc
                        for ww in range(white_list.shape[0]):
                            diff = self.haversine(center_location, white_list[ww])
                            if diff < 5000:
                                ignore = False
                                break
                    if ignore:
                        continue

                # don't add the contour if it's outside of the figure
                indexes = np.where((path.vertices[:,0] > self.S.xlim[0]) &
                                   (path.vertices[:,0] < self.S.xlim[1]) &
                                   (path.vertices[:,1] > self.S.ylim[0]) &
                                   (path.vertices[:,1] < self.S.ylim[1]))
                if len(indexes[0]) == 0:
                    continue

                filtered_paths.append(path)


            print("adding",len(filtered_paths),"filtered paths")

            # important to sort by order first to do containing metrics
            filtered_paths.sort(key = lambda x : len(x.vertices))

            contour_paths.append(filtered_paths)

        return contour_paths


    def remove_empty_df(self, df):
        """Removes empty dataframe objects

        Parameters
        ----------
        df : list
            List of dataframes

        Returns
        -------
        df_new : list
            List of dataframes with empty dataframes removed.

        """
        df_new = []
        old_length = len(df)
        for dd, df in enumerate(df):
            if len(df) > 0:
                df_new.append(df)

        new_length = len(df_new)
        print("removed",old_length-new_length,"empty dataframes")

        return df_new

    def save_and_close(self, figure, plot_name = "fig", type = "png"):
        """Save and close figures.

        Also makes adjustments so figures are displayed neatly. I don't
        know if all of these adjustments are strictly necessary to
        remove the blank space around the figure, but I haven't found
        a different combination that also works. See this post for more
        info on the margins/whitespace adjustments.
        https://stackoverflow.com/a/27227718/12995548

        Parameters
        ----------
        figure : Figure
            Figure object that needs to be saved
        plot_name : string
            Plot name that will be added to file name.

        """

        fig = plt.figure(figure.number)
        ax = plt.subplot(1,1,1)

        # proof fig
        ax.set_xlim(self.S.xlim)
        ax.set_ylim(self.S.ylim)

        # remove axis ticks and labels
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)
        plt.margins(0,0)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())

        fig.savefig(os.path.join(self.repo_dir,"plots",
                    self.S.name + "-" + plot_name + "-" \
                  + self.plot_time + "." + type),
                    format=type,
                    dpi=300.,
                    pad_inches=0,
                    bbox_inches="tight")
        plt.close(fig)


    def laser_contours(self):
        """Laser contours for the entire Bay area.

        """
        level_indexes = [3, 2, 1]

        df_temp = copy.deepcopy(self.df)

        # contour figures in order of [1st level, 2nd level, etc.]
        # cfigs = [ for fi in level_indexes]
        # path figures
        pfigs = [plt.figure(figsize=self.S.figsize) for fi in level_indexes]

        contour_paths = self.get_contour_paths()

        for cc, contour_level in enumerate(contour_paths):
            contour_fig = plt.figure(figsize=self.S.figsize)
            ax = plt.subplot(1,1,1)
            for path in contour_level:
                patch = patches.PathPatch(path, edgecolor='k', facecolor='w', lw=0.072)
                ax.add_patch(patch)
            self.save_and_close(contour_fig,"contour-" + str(cc + 1), "svg")

        print("Done with all contour vectors.\n")

        ################################################################
        # add paths proof fig
        ################################################################

        path_fig = plt.figure(figsize=self.S.figsize)
        path_ax = plt.subplot(1,1,1)

        percent10 = int(len(df_temp)/10.)
        for dd, df in enumerate(df_temp):
            print("dd:",dd,"/",len(df_temp))
            if (dd+1) % percent10 == 0:
                print(math.ceil(100*float(dd)/len(df_temp)),"% of activities plotted")

            if len(df) == 0:
                continue

            if len(df) > self.S.path_min:
                path_ax.plot(df["longitude"].to_numpy().reshape(-1,1),
                         df["latitude"].to_numpy().reshape(-1,1),
                         "k",linewidth = 0.5)

        # proof fig
        path_ax.set_xlim(self.S.xlim)
        path_ax.set_ylim(self.S.ylim)

        # remove axis ticks and labels
        path_ax.axes.xaxis.set_visible(False)
        path_ax.axes.yaxis.set_visible(False)

        path_fig.savefig(os.path.join(self.repo_dir,"plots",
                    self.S.name + "-path-proof-" \
                  + self.plot_time + ".svg"),
                    format="svg",
                    pad_inches=0.0,
                    dpi=300.,
                    bbox_inches="tight")
        plt.close(path_fig)
        print("done with all path proof plotting")

        ################################################################
        # add paths to individual figures
        ################################################################
        background_vertices = np.array([[self.S.xlim[0],self.S.ylim[0]],
                                        [self.S.xlim[0],self.S.ylim[1]],
                                        [self.S.xlim[1],self.S.ylim[1]],
                                        [self.S.xlim[1],self.S.ylim[0]],
                                        ])
        background_path = Path(background_vertices)

        for level_index in range(2,-1,-1):
            for pp, path in enumerate(contour_paths[level_index]):
                print("pp:",pp,"/",len(contour_paths[level_index]))
                level = level_index
                for cp, check_path in enumerate(contour_paths[level_index]):
                    contained_measure = sum(check_path.contains_points(path.vertices))/float(path.vertices.shape[0])
                    if contained_measure > self.S.contained_measure and pp != cp:
                        level = max(0, level_index - 1)
                        print("\n breaking!")
                        # new_fig = plt.figure()
                        # plt.plot(path.vertices[:,0],path.vertices[:,1],"b")
                        # plt.plot(check_path.vertices[:,0],check_path.vertices[:,1],"r")
                        # plt.show()
                        break

                df_cropped = []
                for dd, df in enumerate(df_temp):
                    latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
                                        df["latitude"].to_numpy().reshape(-1,1)))

                    mask = path.contains_points(latlon)
                    mask_idx = np.where(mask)[0]

                    if len(mask_idx) > 0:
                        longs = df.loc[mask_idx,"longitude"].to_numpy().reshape(-1,1)
                        lats = df.loc[mask_idx,"latitude"].to_numpy().reshape(-1,1)
                        latlon2 = np.hstack((lats,longs))
                        separated_latlon = self.separate(latlon2)
                        for sep_latlon in separated_latlon:
                            # fig = plt.figure(cfigs[level].number)
                            # ax = plt.subplot(1,1,1)
                            # ax.plot(sep_latlon[:,1], sep_latlon[:,0], c="b", linewidth = 0.5)
                            fig = plt.figure(pfigs[level].number)
                            ax = plt.subplot(1,1,1)
                            ax.plot(sep_latlon[:,1], sep_latlon[:,0], c="b", linewidth = 0.072)


                    df.drop(labels=np.where(mask)[0], inplace=True)
                    df.reset_index(inplace=True, drop=True)

                    df_cropped.append(df)

                df_temp = self.remove_empty_df(df_cropped)

        # lastly, check whether it's in the water
        for dd, df in enumerate(df_temp):
            level = 0
            latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
                                df["latitude"].to_numpy().reshape(-1,1)))

            mask = background_path.contains_points(latlon)
            mask_idx = np.where(mask)[0]

            if len(mask_idx) > 0:
                longs = df.loc[mask_idx,"longitude"].to_numpy().reshape(-1,1)
                lats = df.loc[mask_idx,"latitude"].to_numpy().reshape(-1,1)
                latlon2 = np.hstack((lats,longs))
                separated_latlon = self.separate(latlon2)
                for sep_latlon in separated_latlon:
                    # fig = plt.figure(cfigs[level].number)
                    # ax = plt.subplot(1,1,1)
                    # ax.plot(sep_latlon[:,1], sep_latlon[:,0], "b", linewidth = 0.5)
                    fig = plt.figure(pfigs[level].number)
                    ax = plt.subplot(1,1,1)
                    ax.plot(sep_latlon[:,1], sep_latlon[:,0], "b", linewidth = 0.072)



        for level_index in level_indexes:
            # save and close all path plots
            fig = plt.figure(pfigs[level_index - 1].number)
            plt.xlim(self.S.xlim)
            plt.ylim(self.S.ylim)

            # remove axis ticks and labels
            ax = plt.gca()
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)

            fig.savefig(os.path.join(self.repo_dir,"plots",
                        self.S.name + "-path-" \
                      + self.plot_time \
                      + "-" + str(level_index) + ".svg"),
                        format="svg",
                        dpi=300.,
                        pad_inches=0.0,
                        bbox_inches="tight",
                        )
            plt.close(pfigs[level_index - 1])

        ################################################################
        # visual contour proof
        ################################################################
        proof_fig = plt.figure(figsize=self.S.figsize)
        proof_ax = plt.subplot(1,1,1)


        background_vertices = np.array([[self.S.xlim[0],self.S.ylim[0]],
                                        [self.S.xlim[0],self.S.ylim[1]],
                                        [self.S.xlim[1],self.S.ylim[1]],
                                        [self.S.xlim[1],self.S.ylim[0]],
                                        ])
        background_path = Path(background_vertices)
        background = patches.PathPatch(background_path, edgecolor='k',
                                        facecolor=self.S.colors[0], lw=2)
        proof_ax.add_patch(background)

        for level_index in range(3):
            """Add new paths."""

            hole_patches = []
            for pp, path in enumerate(contour_paths[level_index]):

                for cp, check_path in enumerate(contour_paths[level_index]):
                    contained_measure = sum(check_path.contains_points(path.vertices))/float(path.vertices.shape[0])
                    # print(check_path.contains_points(path.vertices))
                    # print(contained_measure)
                    if contained_measure > self.S.contained_measure and pp != cp:
                        color = self.S.colors[level_index]
                        #
                        # print(pp,cp,"hey",color,level_index)
                        # ax_new = plt.gca()
                        # # proof fig
                        # ax_new.set_xlim(self.S.xlim)
                        # ax_new.set_ylim(self.S.ylim)
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
                        color = self.S.colors[level_index + 1]

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

            if len(df) > self.S.path_min:
                proof_ax.plot(df["longitude"].to_numpy().reshape(-1,1),
                         df["latitude"].to_numpy().reshape(-1,1),
                         "k",linewidth = 0.5)
            # else:
            #     plt.plot(df["longitude"].to_numpy().reshape(-1,1),
            #     df["latitude"].to_numpy().reshape(-1,1),
            #     "b",linewidth = 0.3)


        # proof fig
        proof_ax.set_xlim(self.S.xlim)
        proof_ax.set_ylim(self.S.ylim)

        # remove axis ticks and labels
        proof_ax.axes.xaxis.set_visible(False)
        proof_ax.axes.yaxis.set_visible(False)
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)
        plt.margins(0,0)
        proof_ax.xaxis.set_major_locator(plt.NullLocator())
        proof_ax.yaxis.set_major_locator(plt.NullLocator())

        proof_fig.savefig(os.path.join(self.repo_dir,"plots",
                    self.S.name + "-contour-proof-" \
                  + self.plot_time + ".png"),
                    format="png",
                    dpi=300.,
                    pad_inches=0,
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

    def separate(self, latlon):
        """Separate out lat lons if they're far apart.

        """
        separated = []
        latlon1 = latlon[:-1,:]
        latlon2 = latlon[1:,:]

        distances = self.haversine_vectorized(latlon1, latlon2)

        disjoint_indexes = np.where(distances > 50.)[0]
        if len(disjoint_indexes) > 0:
            start_index = 0
            for disjoint_idx in disjoint_indexes:
                end_index = disjoint_idx + 1
                separated.append(latlon[start_index:end_index,:])
                start_index = end_index
            if end_index < latlon.shape[0] - 1:
                separated.append(latlon[end_index:,:])
        else:
            separated.append(latlon)

        return separated

    def haversine_vectorized(self, coord1, coord2):
        """
        Taken from: https://janakiev.com/blog/gps-points-distance-python/
        """
        R = 6372800  # Earth radius in meters
        lat1 = coord1[:,0]
        lon1 = coord1[:,1]
        lat2 = coord2[:,0]
        lon2 = coord2[:,1]

        phi1, phi2 = np.radians(lat1), np.radians(lat2)
        dphi       = np.radians(lat2 - lat1)
        dlambda    = np.radians(lon2 - lon1)

        a = np.sin(dphi/2)**2 + \
            np.cos(phi1)*np.cos(phi2)*np.sin(dlambda/2)**2

        return 2*R*np.arctan2(np.sqrt(a), np.sqrt(1 - a))
