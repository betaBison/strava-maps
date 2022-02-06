"""Contour numpy settings

"""

__authors__ = "D. Knowles"
__date__ = "31 Jan 2022"

import time
import numpy as np

class Tristan():
    """Tristan's preferences. """

    def __init__(self):
        self.name = "tristan"
        self.levels = [-1000, 1.5, 100., 250, 5000]
        self.colors = ["#BDD1E8","#6c5a47","#865c3c","#c08e5f"]
        self.xlim = [-122.78779982641429,-121.67751846782933]
        self.ylim = [36.94154503938612, 38.04094280129413]
        # see full page
        # self.xlim = [-123.78779982641429,-120.67751846782933]
        # self.ylim = [35.94154503938612, 39.04094280129413]
        self.figsize = (15.75, 19.75)
        # self.plot_name = "tristan-" + str(time.time())[:10] + ".png"
        self.blur = 1
        self.contour_min = 110
        self.path_min = 400
        # white list contours
        self.wlc = np.array([[37.86171298036722, -122.43017188317853], # angel island
                             [37.82634816468346, -122.4226103787972], # alcatrz island
                             [37.82438143574645, -122.37045086434782], # treasure island
                             [37.810816442739046, -122.36563221855539], # yerba buena island
                             ])
        self.contained_measure = 0.9

class Derek():
    """Derek's preferences. """

    def __init__(self):
        self.name = "derek"
        self.levels = [-1000, 3, 150., 350., 5000]
        self.colors = ["#BDD1E8","#6c5a47","#865c3c","#c08e5f"]
        self.xlim = [-122.78779982641429,-121.67751846782933]
        self.ylim = [36.94154503938612, 38.04094280129413]
        self.figsize = (15.75, 19.75)
        self.blur = 0.1
        self.contour_min = 110
        self.path_min = 400
        # white list contours
        self.wlc = np.array([[37.86171298036722, -122.43017188317853], # angel island
                             [37.82634816468346, -122.4226103787972], # alcatrz island
                             [37.82438143574645, -122.37045086434782], # treasure island
                             [37.810816442739046, -122.36563221855539], # yerba buena island
                             ])
        self.contained_measure = 0.9
