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
from tqdm import tqdm
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
        activities_path = os.path.join(data_path,"activities_csv")

        file_dir = os.path.dirname(os.path.realpath(__file__))
        self.repo_dir = os.path.join(file_dir, "..")

        self.num_levels = len(self.S.levels) - 2 # number of levels

        self.df = [] # holds all activity data

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

        # create background path for plotting
        background_vertices = np.array([[self.S.xlim[0],self.S.ylim[0]],
                                        [self.S.xlim[0],self.S.ylim[1]],
                                        [self.S.xlim[1],self.S.ylim[1]],
                                        [self.S.xlim[1],self.S.ylim[0]],
                                        ])
        self.background_path = Path(background_vertices)

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
        if self.S.debug:
            self.save_and_close(temp_fig, "debug-raw-contours")
        else:
            plt.close(temp_fig)

        return contours

    def get_contour_paths(self):
        """Get contour paths for all levels.

        Returns the collection of contours as paths.

        Also removes contour points outside of bounds and return as a
        collection of numpy arrays.

        Returns
        -------
        contour_paths : list
            List of lists of all contour paths.
            Order of [[1st level paths],[2nd level paths],[etc...]].
        contour_vertices : list
            List of lists of numpy arrays representing contour vertices.

        """

        contour_paths = []
        contour_vertices = []
        contours = self.get_contours()

        for level_index in range(1, self.num_levels + 1):
            # remove small collections unless they're in the white list
            filtered_paths = []
            filtered_vertices = []
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
                                   (path.vertices[:,1] < self.S.ylim[1]))[0]

                if len(indexes) != 0:
                    # add Path objects
                    filtered_paths.append(path)
                    # also add vertices objects
                    # separate out into distict line segments if far apart
                    separated_latlon = self.separate(np.flip(path.vertices[indexes],axis=1),
                                                     self.S.contour_separate)
                    for sep_latlon in separated_latlon:
                        # flip back to be long, lat as normal with vertices
                        filtered_vertices.append(np.flip(sep_latlon, axis=1))


            print("adding",len(filtered_paths),"filtered paths")

            # important to sort by order first to do containing metrics
            filtered_paths.sort(key = lambda x : len(x.vertices))

            contour_paths.append(filtered_paths)
            contour_vertices.append(filtered_vertices)

        return contour_paths, contour_vertices

    def check_level_containment(self, raw_contour_paths):
        """Restructures contours to their respective level.

        If a contour path is within another contour path of it's same
        levle, than it is a valley and is moved down a level.

        Parameters
        -------
        raw_contour_paths : list
            List of lists of all contour paths.
            Order of [[1st level paths],[2nd level paths],[etc...]].

        Returns
        -------
        check_contours : list
            List of lists of all contour paths but grouped by level.
            Order of [[1st level paths],[2nd level paths],[etc...]].

        """
        check_contours = [[] for ll in range(self.num_levels*2)]
        for level_index in range(self.num_levels):
            for pp, path in enumerate(raw_contour_paths[level_index]):
                level = 2*level_index
                for cp, check_path in enumerate(raw_contour_paths[level_index]):
                    contained_measure = sum(check_path.contains_points(path.vertices))/float(path.vertices.shape[0])
                    if contained_measure > self.S.contained_measure and pp != cp:
                        level = 2*level_index + 1
                        break
                check_contours[level].append(path)

        return check_contours


    def get_activity_vertices(self, check_contours):
        """Returns all activity vertices in level sets.

        Parameters
        ----------
        check_contours : list
            List of lists of all contour paths but grouped by level.
            Order of [[1st level paths],[2nd level paths],[etc...]].

        Returns
        -------
        contour_vertices : list
            List of lists of numpy arrays representing contour vertices.

        """
        activity_vertices = [[] for ll in range(self.num_levels)]

        # path figures
        df_temp = copy.deepcopy(self.df)

        for level in tqdm(range(len(check_contours)-1,-1,-1),desc="contour levels:",position=0):
            for pp, path in enumerate(tqdm(check_contours[level],desc="paths",position=1,leave=False)):
                if level % 2 == 1:
                    plot_level = int((level - 3)/2)
                else:
                    plot_level = int(level/2)
                plot_level = max(0,plot_level)

                df_cropped = []
                for dd, df in enumerate(tqdm(df_temp,desc="activities",position=2,leave=False)):
                    latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
                                        df["latitude"].to_numpy().reshape(-1,1)))

                    mask = path.contains_points(latlon)
                    # also check that activity is within bounds
                    back_mask = self.background_path.contains_points(latlon)
                    mask_idx = np.where(mask & back_mask)[0]
                    # mask_idx = np.where(mask)[0]

                    if len(mask_idx) > 0:
                        longs = df.loc[mask_idx,"longitude"].to_numpy().reshape(-1,1)
                        lats = df.loc[mask_idx,"latitude"].to_numpy().reshape(-1,1)
                        latlon = np.hstack((lats,longs))
                        if len(activity_vertices[plot_level]) == 0:
                            activity_vertices[plot_level].append(latlon)
                        else:
                            num_new = latlon.shape[0]
                            threshold_vals = np.zeros((num_new,))
                            for previous_points in activity_vertices[plot_level]:
                                num_prev = previous_points.shape[0]
                                prev_tile = np.tile(previous_points,(num_new,1))
                                new_tile = np.repeat(latlon,num_prev,axis=0)

                                new_dist = self.haversine_vectorized(prev_tile,new_tile).reshape(num_new,num_prev)
                                sub_vals = np.where(new_dist > self.S.activity_redundancy, 0, 1)
                                sub_vals = np.sum(sub_vals,axis=1)
                                threshold_vals += sub_vals

                            to_add = latlon[np.where(threshold_vals < 1)[0],:]
                            activity_vertices[plot_level].append(to_add)

                    df.drop(labels=np.where(mask)[0], inplace=True)
                    df.reset_index(inplace=True, drop=True)

                    df_cropped.append(df)

                df_temp = self.remove_empty_df(df_cropped)

        # lastly, check whether it's in the water
        for dd, df in enumerate(df_temp):
            plot_level = 0
            latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
                                df["latitude"].to_numpy().reshape(-1,1)))

            mask = self.background_path.contains_points(latlon)
            mask_idx = np.where(mask)[0]

            if len(mask_idx) > 0:
                longs = df.loc[mask_idx,"longitude"].to_numpy().reshape(-1,1)
                lats = df.loc[mask_idx,"latitude"].to_numpy().reshape(-1,1)
                latlon = np.hstack((lats,longs))
                if len(activity_vertices[plot_level]) == 0:
                    activity_vertices[plot_level] = latlon
                else:
                    num_new = latlon.shape[0]
                    threshold_vals = np.zeros((num_new,))
                    for previous_points in activity_vertices[plot_level]:
                        num_prev = previous_points.shape[0]
                        prev_tile = np.tile(previous_points,(num_new,1))
                        new_tile = np.repeat(latlon,num_prev,axis=0)

                        new_dist = self.haversine_vectorized(prev_tile,new_tile).reshape(num_new,num_prev)
                        sub_vals = np.where(new_dist > self.S.activity_redundancy, 0, 1)
                        sub_vals = np.sum(sub_vals,axis=1)
                        threshold_vals += sub_vals

                    to_add = latlon[np.where(threshold_vals < 1)[0],:]
                    activity_vertices[plot_level].append(to_add)

        separated_vertices = []
        for ll in range(self.num_levels):
            separated_vertices.append([])
            for activity in activity_vertices[ll]:
                separated_latlon = self.separate(activity,
                                    self.S.activity_separate)
                for sep_latlon in separated_latlon:
                    if len(sep_latlon) > 20:
                        separated_vertices[ll].append(np.flip(sep_latlon,axis=1))

        return separated_vertices


    def remove_empty_df(self, df, verbose = False):
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
        if verbose:
            print("removed",old_length-new_length,"empty dataframes")

        return df_new


    def crop_df(self, df, verbose = False):
        """Removes dataframe objects that are outside the fig bounds.

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
            latlon = np.hstack((df["longitude"].to_numpy().reshape(-1,1),
                                df["latitude"].to_numpy().reshape(-1,1)))
            back_mask = self.background_path.contains_points(latlon)
            mask_idx = np.where(back_mask)[0]
            # mask_idx = np.where(mask)[0]

            if len(mask_idx) > 0:
                df_new.append(df)

        new_length = len(df_new)
        if verbose:
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


        self.df = self.crop_df(self.df)

        contour_paths, contour_vertices = self.get_contour_paths()

        print("plotting contour vertices")
        for cv, contour_level in enumerate(contour_vertices):
            contour_fig = plt.figure(figsize=self.S.figsize)
            ax = plt.subplot(1,1,1)
            if self.S.debug:
                contour_debug = plt.figure(figsize=self.S.figsize)
                debug_ax = plt.subplot(1,1,1)

            for vertices in contour_level:
                ax.plot(vertices[:,0], vertices[:,1], c="k",
                        linewidth = self.S.laser_linewidth)
                if self.S.debug:
                    # allows you to tell where contours are separated
                    debug_ax.plot(vertices[:,0], vertices[:,1],
                                  linewidth = 3.)
            self.save_and_close(contour_fig,"contour-" + str(cv + 1), "svg")
            if self.S.debug:
                self.save_and_close(contour_debug,
                        "debug-contour-vertices-" + str(cv + 1), "svg")

        print("Done with all contour vectors.\n")

        ################################################################
        # add activity paths to individual figures
        ################################################################

        print("checking level containment")
        check_contours = self.check_level_containment(contour_paths)
        print("getting activity vertices")
        activity_vertices = self.get_activity_vertices(check_contours)

        print("plotting activty vertices")

        for al, activity_level in enumerate(activity_vertices):
            activity_fig = plt.figure(figsize=self.S.figsize)
            ax = plt.subplot(1,1,1)
            if self.S.debug:
                activity_debug = plt.figure(figsize=self.S.figsize)
                debug_ax = plt.subplot(1,1,1)

            for vertices in activity_level:
                ax.plot(vertices[:,0], vertices[:,1], c="k",
                        linewidth = self.S.laser_linewidth)
                if self.S.debug:
                    # allows you to tell where contours are separated
                    debug_ax.plot(vertices[:,0], vertices[:,1],
                                  linewidth = 3.)
            self.save_and_close(activity_fig,"activity-" + str(al + 1), "svg")
            if self.S.debug:
                self.save_and_close(activity_debug,
                        "debug-activity-vertices-" + str(al + 1), "svg")

        print("Done with all activity path plots.")

        ################################################################
        # visual contour proof
        ################################################################
        proof_fig = plt.figure(figsize=self.S.figsize)
        proof_ax = plt.subplot(1,1,1)

        background = patches.PathPatch(self.background_path, edgecolor='k',
                                        facecolor=self.S.colors[0], lw=2)

        if self.S.debug:
            proof_debug = plt.figure(figsize=self.S.figsize)
            debug_ax = plt.subplot(1,1,1)
            debug_ax.add_patch(copy.copy(background))

        proof_ax.add_patch(background)



        for level_index in range(self.num_levels):
            """Add new paths."""

            hole_patches = []
            for pp, path in enumerate(contour_paths[level_index]):

                for cp, check_path in enumerate(contour_paths[level_index]):
                    contained_measure = sum(check_path.contains_points(path.vertices))/float(path.vertices.shape[0])
                    if contained_measure > self.S.contained_measure and pp != cp:
                        color = self.S.colors[level_index]

                        new_patch = patches.PathPatch(path, edgecolor="k",
                                                        facecolor=color, lw=0)
                        hole_patches.append(new_patch)

                    else:
                        color = self.S.colors[level_index + 1]

                    proof_patch = patches.PathPatch(path, edgecolor="k",
                                                    facecolor = color, lw=0)
                    if self.S.debug:
                        debug_ax.add_patch(copy.copy(proof_patch))
                    proof_ax.add_patch(proof_patch)

            # add hole patches on top
            for hole_patch in hole_patches:
                if self.S.debug:
                    debug_ax.add_patch(copy.copy(hole_patch))
                proof_ax.add_patch(hole_patch)

            print("added proof contours for level",level_index+1,"/3")

        for al, activity_level in enumerate(activity_vertices):
            for vertices in activity_level:
                proof_ax.plot(vertices[:,0], vertices[:,1], c="k",
                        linewidth = 1.0)
                if self.S.debug:
                    debug_ax.plot(vertices[:,0], vertices[:,1], c="C"+str(al),
                            linewidth = 1.0)

        self.save_and_close(proof_fig,
                        "contour-proof", "png")
        if self.S.debug:
            self.save_and_close(proof_debug,
                            "debug-contour-proof", "png")


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

    def separate(self, latlon, threshold):
        """Separate out lat lons if they're far apart.

        Assume latitude is in leftmost column and longitude is in
        rightmost column.

        Parameters
        ----------
        latlon : np.ndarray
            Array of shape [n x 2] where latitude is leftmost column and
            longitude is rightmost column.
        threshold : float
            Distance at which arrays should be separated into distinct
            groups.

        Returns
        -------
        separated : list
            List of separated latlon arrays same as the input.

        """
        separated = []
        latlon1 = latlon[:-1,:]
        latlon2 = latlon[1:,:]

        distances = self.haversine_vectorized(latlon1, latlon2)

        disjoint_indexes = np.where(distances > threshold)[0]
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
        """Calculates distance between points.

        Parameters
        ----------
        coord1 : np.ndarray
            Array of shape [n x 2] where latitude is leftmost column and
            longitude is rightmost column.
        coord2 : np.ndarray
            Array of shape [n x 2] where latitude is leftmost column and
            longitude is rightmost column.

        Returns
        -------
        distance : np.ndarray
            Distances between points as an [n,] size array

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
