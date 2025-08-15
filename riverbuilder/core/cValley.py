'''This is the valley class module.

It handles simulating a valley containing a river.
It consists of a river, which is a Channel object, and arbitrary
number of valley levels.

All offsets from valley levels to river are calculated in xy-coordinate system.
'''

import numpy as np
from . import functions
from .cPipe import Pipe
from .cChannel import Channel
from math import pi, sqrt, log
import csv
import matplotlib.pyplot as plt


class Valley(Pipe):

    def __init__(self, channel=Channel()):
        '''Initiator for Valley class

        x_len -- length of valley in x direction
        channel -- Channel object passed to a valley
        x_slope -- slope in x direction
        dx -- resolution in x direction

        Class attributes:
        reshape -- flag show that if the channel has been reshaped
        channelCenter_x -- x values for the points of channel centerline
        channelCenter_y -- y values for the points of channel centerline
        '''
        self.channel = channel
        super().__init__(int(channel.x_len), channel.x_slope, channel.dx, channel.zd)
        self.channelCenter_x = channel.x_v
        self.channelCenter_y = channel.getCenterline_y()
        self.setThalweg()


    # def setValleyBoundary(self, z_offset, y_offset, direction, yfun=None):
    #     try:
    #         print(f"Called setValleyBoundary for direction: {direction}")
    #         z_start = self.getThalweg() - self.zd
    #         z_max = np.amax(self.channel.levels_z[direction][-1])*self.dx
    #         print(f"----------levels_z[{direction}] = {self.channel.levels_z[direction]}---------")
    #         z_offset += z_max
    #     except Exception as e:
    #         print(f"Error in setValleyBoundary: {e}")

    #     if direction == 'left':
    #         y_max = np.amax(self.channel.levels_y[direction][-1])*self.dx
    #         y_offset += y_max
    #     else:
    #         y_max = np.amin(self.channel.levels_y[direction][-1])*self.dx
    #         y_offset = abs(y_max - y_offset)


    #     super().setLevel(z_offset, z_start, y_offset, direction, yfun)


    def setValleyBoundary(self, z_offset, y_offset, direction, yfun=None):
        # Compute thalweg start
        z_start = self.getThalweg() - self.zd

        #  bank levels:
        if not self.channel.levels_z[direction]:
            raise RuntimeError(f"No elevation levels built for channel.{direction!r}")
        if not self.channel.levels_y[direction]:
            raise RuntimeError(f"No lateral levels built for channel.{direction!r}")

        # Getting  the last bank‐level arrays
        z_arr = self.channel.levels_z[direction][0]
        y_arr = self.channel.levels_y[direction][0]

        # Check they're non‐empty
        if z_arr.size == 0:
            raise RuntimeError(f"channel.levels_z[{direction}][-1] is empty")
        if y_arr.size == 0:
            raise RuntimeError(f"channel.levels_y[{direction}][-1] is empty")

        # Compute max/min offsets
        z_offset += np.amax(z_arr) * self.dx
        if direction == 'left':
            y_offset += np.amax(y_arr) * self.dx
        else:
            y_offset = abs(np.amin(y_arr) * self.dx - y_offset)


        super().setLevel(z_offset, z_start, y_offset, direction, yfun)


    def getLevel(self, dictkey, direction, ind=None):
        '''Rewrite the getLevel function. Use the getBound for the first time.'''
        if self.leveldict[dictkey][direction] == []:
            return self.channel.getLevel(dictkey, direction, ind)
        else:
            return super().getLevel(dictkey, direction, ind)


    def tolist_levelxy(self):
        '''Return x, y values for all levels in a secondary list.'''
        x = [] + self.channelCenter_x.tolist()
        y = [] + self.channelCenter_y.tolist()

        self.helpAppend(x, self.channel.levels_x)
        self.helpAppend(y, self.channel.levels_y)
        self.helpAppend(x, self.levels_x)
        self.helpAppend(y, self.levels_y)
        return [x, y]

    def tolabeledlist_levelxy(self):
        '''Return x, y values for all levels in a secondary list.'''
        x = []
        y = []
        label = []

        # Add left valley
        for i in range(len(self.levels_x['left'])):
            x += self.levels_x['left'][i].tolist()
            y += self.levels_y['left'][i].tolist()
            label += ['LV'+str(i)]*len(self.levels_x['left'][i])

        # Add left channel bank
        for i in range(len(self.channel.levels_x['left'])):
            x += self.channel.levels_x['left'][i].tolist()
            y += self.channel.levels_y['left'][i].tolist()
            label += ['LB'+str(i)]*len(self.channel.levels_x['left'][i])

        # Add channel center line
        x += self.channelCenter_x.tolist()
        y += self.channelCenter_y.tolist()
        label += ['CL'] * len(self.channelCenter_x.tolist())

        # Add right channel bank
        for i in range(len(self.channel.levels_x['right'])):
            x += self.channel.levels_x['right'][i].tolist()
            y += self.channel.levels_y['right'][i].tolist()
            label += ['RB'+str(i)]*len(self.channel.levels_x['right'][i])

        # Add right valley
        for i in range(len(self.levels_x['right'])):
            x += self.levels_x['right'][i].tolist()
            y += self.levels_y['right'][i].tolist()
            label += ['RV'+str(i)]*len(self.levels_x['right'][i])

        return [x, y, label]

    def tolist_levelxz(self):
        '''Return x, y values for all levels in a secondary list.'''
        x = [] + self.channelCenter_x.tolist()
        z = [] + self.channel.getThalweg().tolist()

        self.helpAppend(x, self.channel.levels_x)
        self.helpAppend(z, self.channel.levels_z)
        self.helpAppend(x, self.levels_x)
        self.helpAppend(z, self.levels_z)

        return [x, z]


    def tolist_all(self):
        '''Return x, y z values for everythin in a secondary list.'''

        xshape_x, xshape_y, xshape_z = self.channel.getXShape()
        center_z = self.channel.getCenterlineElevation()
        x = [] + self.channelCenter_x.tolist()
        y = [] + self.channelCenter_y.tolist()
        z = [] + center_z.tolist()
        label = ['CL'] * len(x)

        x += xshape_x.tolist()
        y += xshape_y.tolist()
        z += xshape_z.tolist()
        label += ['XS']*len(xshape_x)

        self.helpAppend(x, self.channel.levels_x)
        self.helpAppend(y, self.channel.levels_y)
        self.helpAppend(z, self.channel.levels_z)
        label += ['RB'] * (len(x) - len(label))

        self.helpAppend(x, self.levels_x)
        self.helpAppend(y, self.levels_y)
        self.helpAppend(z, self.levels_z)
        label += ['VL'] * (len(x) - len(label))


        return [x, y, z, label]


    def getXShapePlot(self):
        """
        return matplotlib plot object that contains X-Shape plots of the valley,
        merged with the channel cross-section.
        """

        # Finding station (minInd) where slope is min
        slope = self.getSlope()
        minInd = np.argmin(np.abs(slope))
        minInd_x = self.x_v[minInd]
        
        y = []
        z = []

        #Valley LEFT side (in reverse order) – nearest-x lookup
        for i in range(len(self.levels_y['left'])):
            lvl_x = self.levels_x['left'][-1 - i]
            if lvl_x.size == 0:
                continue
            lvl_y = self.levels_y['left'][-1 - i]
            lvl_z = self.levels_z['left'][-1 - i]
            ind_local = np.argmin(np.abs(lvl_x - minInd_x))
            y.append(lvl_y[ind_local] * self.dx)
            z.append(lvl_z[ind_local] * self.dx)

        #Channel LEFT side -- also in reverse order
        for i in range(len(self.channel.levels_y['left'])):
            lvl_x = self.channel.levels_x['left'][-1 - i]
            if lvl_x.size == 0:
                continue
            ind_local = np.argmin(np.abs(lvl_x - minInd_x))
            y.append(self.channel.levels_y['left'][-1 - i][ind_local] * self.dx)
            z.append(self.channel.levels_z['left'][-1 - i][ind_local] * self.dx)

        # find wbf + station indices
        indLeft  = np.argmin(np.abs(self.channel.levels_x['left'][0] - minInd_x))
        indRight = np.argmin(np.abs(self.channel.levels_x['right'][0] - minInd_x))
        wbf = (self.channel.levels_y['left'][0][indLeft]
            - self.channel.levels_y['right'][0][indRight])
        ind = indLeft

        # (PY, CF, AF, DT, TU, or AU/SU)
        ctype = getattr(self.channel, 'cross_section_type', None)
        if ctype == 'CF':
            CF_a = getattr(self.channel, 'CF_a', 5)
            CF_b = getattr(self.channel, 'CF_b', 5)
            CF_c = getattr(self.channel, 'CF_c', 5)
            channel_y, channel_z = self.channel.cfXShape(
                wbf,
                n=self.channel.xshapePoints,
                CF_a=CF_a, CF_b=CF_b, CF_c=CF_c
            )
            channel_y = channel_y[::-1]
            channel_z = channel_z[::-1]

        elif ctype == 'PY':
            channel_y, channel_z = self.channel.pyXShape(wbf)
            channel_y = channel_y[::-1]
            channel_z = channel_z[::-1]

        elif ctype == 'AF':
            d1   = getattr(self.channel, 'af_d1', 5)
            d2   = getattr(self.channel, 'af_d2', 5)
            ang1 = getattr(self.channel, 'af_ang1', 55)
            ang2 = getattr(self.channel, 'af_ang2', 55)
            channel_y, channel_z = self.channel.afXShape(
                wbf,
                n=self.channel.xshapePoints,
                d1=d1, d2=d2, ang1=ang1, ang2=ang2
            )
            channel_y = channel_y[::-1]
            channel_z = channel_z[::-1]

            actual_left_raw = self.channel.levels_z['left'][0][indLeft]
            delta = actual_left_raw - channel_z[0]
            channel_z = channel_z + delta

            centerOffset = (
                self.channel.levels_y['left'][0][indLeft]
            + self.channel.levels_y['right'][0][indRight]
            ) / 2.0
            channel_y = (channel_y + centerOffset) * self.dx
            channel_z = channel_z * self.dx

            y += channel_y.tolist()
            z += channel_z.tolist()
            y.append(np.nan);  z.append(np.nan)

        elif ctype == 'DT':
            channel_y, channel_z = self.channel.dtXShape(wbf)
            channel_y = channel_y[::-1]
            channel_z = channel_z[::-1]

        elif ctype == 'TU':
            channel_y, channel_z = self.channel.tuXShape(wbf)
            channel_y = channel_y[::-1]
            channel_z = channel_z[::-1]

        else:
            if self.channel.tz == -1:
                channel_y, channel_z = self.channel.pointXShape(
                    ind, 0, wbf, self.channel.xshapePoints
                )
            else:
                channel_y, channel_z = self.channel.suXShape(
                    ind, wbf, self.channel.tz, self.channel.xshapePoints
                )

        # Shift channel Y by center offset, then scale
        centerOffset = (
            self.channel.levels_y['left'][0][indLeft]
        + self.channel.levels_y['right'][0][indRight]
        ) / 2.0
        channel_y = (channel_y + centerOffset) * self.dx
        channel_z = channel_z * self.dx

        # Append the channel cross-section
        y += channel_y.tolist()
        z += channel_z.tolist()

        # Channel RIGHT side (forward order)
        for i in range(len(self.channel.levels_y['right'])):
            lvl_x = self.channel.levels_x['right'][i]
            if lvl_x.size == 0:
                continue
            ind_local = np.argmin(np.abs(lvl_x - minInd_x))
            y.append(self.channel.levels_y['right'][i][ind_local] * self.dx)
            z.append(self.channel.levels_z['right'][i][ind_local] * self.dx)

        # Valley RIGHT side (forward order) – nearest-x lookup
        for i in range(len(self.levels_y['right'])):
            lvl_x = self.levels_x['right'][i]
            if lvl_x.size == 0:
                continue
            lvl_y = self.levels_y['right'][i]
            lvl_z = self.levels_z['right'][i]
            ind_local = np.argmin(np.abs(lvl_x - minInd_x))
            y.append(lvl_y[ind_local] * self.dx)
            z.append(lvl_z[ind_local] * self.dx)

        # Plot
        fig, ax = plt.subplots(1, 1)
        fig.suptitle('Valley X-Shape')
        ax.plot(y, z, '-', marker='o')
        ax.set_xlabel('Y')
        ax.set_ylabel('Z')
        return fig




    def tocsv(self, outfile):
        '''Outwrite xy values and xz values of all levels to output file.'''
        header = ["X", "Y", "Z", "Label"]
        out = [header]
        xyz = self.tolist_all()

        xyz_out = [[round(xyz[0][i]*self.dx,3), round(xyz[1][i]*self.dx, 3), round(xyz[2][i]*self.dx, 3), xyz[3][i]] for i in range(len(xyz[0]))]
#        xyz_out = [[round(xyz[0][i],3), round(xyz[1][i], 3), round(xyz[2][i], 3), xyz[3][i]] for i in range(len(xyz[0]))]
        out += xyz_out

        with open(outfile+".csv", 'w') as cf:
            cw = csv.writer(cf, lineterminator='\n')
            cw.writerows(out)


    def getValleySlope(self):
        '''Return Valley slope.'''
        return self.getPipeSlope()


    def __str__(self):
        '''Turn most important data of valley.'''
        sl = self.getSL()
        slope = self.getValleySlope()
        s = 'Sinuosity:'+str(round(sl, 3))+'\n'
        s += 'Valley Slope:'+str(slope)+'\n'

        for i in range(0, len(self.levels_n['left'])):
            s += 'Average Width of '+'L'+str(i)+' Valley Level is: '+str(round(np.average(self.levels_n['left'][i])*self.dx,3)) + '\n'
            
        for i in range(0, len(self.levels_n['right'])):
            s += 'Average Width of '+'R'+str(i)+' Valley Level is: '+str(abs(round(np.average(self.levels_n['right'][i]*self.dx), 3))) + '\n'

        return s
###########################################################

    def helpAppend(self, li, dic):
        for array in dic["left"]:
            li += array.tolist()
        for array in dic["right"]:
            li += array.tolist()