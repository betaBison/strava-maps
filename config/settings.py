"""User settings.

"""

__authors__ = "D. Knowles"
__date__ = "31 Jan 2022"

import time
import numpy as np

class Template():
    """Template preferences. """

    def __init__(self):
        self.name = "your_name"
        # absolute path of the export_XXXXXXXX/ directory you received by
        # exporting all of your Strava data in bulk.
        self.data_path = "/home/user/data/strava/export_12345678/"
        self.topo_name = "SF_bay_topo_medium.npy"
        self.topo_xlim = (-124.000555556493, -120.99944444380549)
        self.topo_ylim = (35.999444443906840, 39.0005555565946)
        self.levels = [-np.inf, 0., 100., 200, np.inf]
        self.colors = ["#BDD1E8","#6c5a47","#865c3c","#c08e5f"]
        self.xlim = [-122.78779982641429,-121.67751846782933]
        self.ylim = [36.94154503938612, 38.04094280129413]
        self.figsize = (16, 20)
        self.blur = 1
        self.contour_min = 200
        self.path_min = -np.inf
        # white list contours
        self.wlc = np.array([[37.86171298036722, -122.43017188317853], # angel island
                             [37.82634816468346, -122.4226103787972], # alcatrz island
                             [37.82438143574645, -122.37045086434782], # treasure island
                             [37.810816442739046, -122.36563221855539], # yerba buena island
                             ])
        self.contained_measure = 0.9
