'''This is the Channel module

It can simulate a river channel basing on the inputs it is provided.
It consists of a centerline, an inner channel, and arbitrary number of
outer banks.
All functions apply to it should be continuous.
The offsets from banks to centerline are in sn coordinate system, and 
transform into xy coordinate system later on.
'''

import numpy as np
from . import functions
import random
import math
from math import pi, sqrt, log, ceil, floor
import csv
from .cPipe import Pipe
import matplotlib.pyplot as plt
import sys

class Channel(Pipe):

    def __init__(self, x_len=100, wbf_min=0, valley_slope=0.01, dx=1, zd=1000, cross_section_type ="default"):
        '''Channel class initiator

        x_len -- int; valley length in x direction
        wbf_min -- float; minimum bankfull width
        valley_slope -- float; slope of valley
        dx -- int; resolution in x direction

        class private variables:
        hbf -- float; average bankfull height
        thalweg -- array; z values of thalweg
        curvature -- array; curvature of centerline
        xshapePoints -- int; number of points in each Xshape
        xshape_x -- array; x values of Xshape
        xshape_y -- array; y values of Xshape
        xshape_z -- array; z values of Xshape
        z_center -- array; z values of centerline
        dynamicCurv -- array; values of curvature of center line
        tz -- int; trapezoid xshape bottom points. -1 means asymetric
        '''
        super().__init__(int(x_len), valley_slope, dx, zd)
        self.wbf_min = wbf_min/dx
        self.cross_section_type = cross_section_type
        self.turns_center = []
        self.hbf = None
        self.curvature = None
        self.xshapePoints = 21
        self.xshape_x = None
        self.xshape_y = None
        self.xshape_z = None
        self.z_center = None
        self.dynamicCurv = None
        self.channelUndulation = None
        self.tz = -1

        self.d1 = 5
        self.d2 = 5
        self.ang1 = 55
        self.ang2 = 55


    def shapeCenterline(self, fun):
        '''Shape the centerline. Basically recalculate the centerline.'''
        x_v = self.x_v
        n_v = self.getCenterline_y()

        x_v_valley = list(set(x_v.tolist()))
        x_v_valley.sort()
        x_v_valley = np.array(x_v_valley)
        x_v_valley, y_v = fun(x_v_valley)
        y_v = y_v/self.dx
        x_max = np.amax(x_v_valley)

        out_x, out_y = [], []
        for i in range(len(x_v)):
            x = x_v[i]
            ind = np.where(x_v_valley == x)[0][0]
            if x == x_max:
                continue
            x1, x2 = x_v_valley[ind], x_v_valley[ind+1]
            y1, y2 = y_v[ind], y_v[ind+1]
            x_new, y_new = functions.sn_to_xy(x1, y1, x2, y2, n_v[i])
            out_x.append(x_new)
            out_y.append(y_new)
        self.x_v = np.array(out_x)
        self.y_center = np.array(out_y)


    def getRiverSlope(self):
        '''Return river slope'''
        return self.getPipeSlope()


    def setXShapePoints(self, n):
        '''Set how many points in one x-section shape.'''
        self.xshapePoints = n


    def setHbfManual(self, hbf):
        '''Mannually set self.hbf'''
        self.hbf = hbf/self.dx


    def setHbf(self, d50=0.01, css=0.047, g_s=1922, g_w=1000):
        '''Automatically calculate Hbf'''
        self.hbf = functions.shields(d50, css, self.getRiverSlope(), g_s, g_w)/self.dx


    def getHbf(self):
        '''Return self.hbf'''
        if self.hbf is None:
            self.setHbf()
        return self.hbf


    def setTZ(self, n):
        self.tz = n


    def setCurvature(self, fun):
        '''Mannually set centerline curvature

        fun -- function to calculate curvature
        '''
        x = self.x_v
        dummy, self.dynamicCurv = fun(x)


    def getDynamicCurv(self):
        '''Return self.dynamicCur'''
        if self.dynamicCurv is None:
            self.setDynamicCurv()
        return self.dynamicCurv


    def createInnerChannel(self, leftFun=None, rightFun=None, thalwegFun=None):
        '''Create most inner channel of river

        leftFun -- function to calculate left inner bank
        rightFun -- function to calculate right innert bank
        thalwegFun -- function to calculate thalweg

        Value will be modified:
        self.levels_x
        self.levels_y
        self.levels_z
        self.levels_n
        '''
        self.setThalweg(thalwegFun)

        thalweg = self.getThalweg()
        orig_thalweg = thalweg+self.channelUndulation
        thalweg_max = np.amax(orig_thalweg)
        z_start = thalweg_max - self.channelUndulation

        hbf = self.getHbf()*self.dx
        wbf = self.wbf_min/2*self.dx

        self.setLevel(hbf, z_start, wbf, 'left', leftFun, True)
        self.setLevel(hbf, z_start, wbf, 'right', rightFun, True)


    def getAveWbf(self):
        '''Return average bankfull width.'''
        if self.levels_y['left'] == []:
            self.createInnerChannel()
        bf = self.levels_n['left'][0] + np.absolute(self.levels_n['right'][0])
        return np.average(bf)*self.dx


    def getAveHbf(self):
        '''Return average bankfull height.'''
        if self.levels_y['left'] == []:
            self.createInnerChannel()
        thalweg = self.getThalweg()
        flat_thalweg = thalweg + self.channelUndulation
        thalweg_max = np.amax(flat_thalweg)
        diff = thalweg_max - flat_thalweg
        return (np.average(diff) + self.getHbf())*self.dx


    def getCoWbf(self):
        '''Return coefficient of variation of bankfull width.'''
        ave = self.getAveWbf()
        std = np.std(self.levels_n['left'][0]*self.dx + (np.absolute(self.levels_n['right'][0])*self.dx))

        return std/ave


    def getCoHbf(self):
        '''Return coefficient of variation of bankfull width.'''
        thalweg = self.getThalweg()
        flat_thalweg = thalweg + self.channelUndulation
        thalweg_max = np.amax(flat_thalweg)
        diff = thalweg_max - flat_thalweg
        ave = (np.average(diff) + self.getHbf())*self.dx
        std = np.std(diff*self.dx)
        return std/ave


    def getXShape(self):
        '''Return x, y, z values for Xshape of the whole channel'''
        if self.xshape_x is None:
            self.setXShape()
        return self.xshape_x, self.xshape_y, self.xshape_z

    def setXShapeAF(self):
        """
        Build the trapezoidal cross-section lines using self.d1, self.d2, etc.
        and store them in self.xshape_x, etc.
        """
        import numpy as np

        out_x, out_y, out_z = [], [], []
        center_z = []
        y_center = self.getCenterline_y()
        s_center = self.getCenterline_sn()
        pipe_slope = self.getPipeSlope()

        xshape_lines = [[] for _ in range(self.xshapePoints)]

        for ind in range(len(y_center) - 1):
            wbf = abs(self.levels_n['left'][0][ind]) + abs(self.levels_n['right'][0][ind])
            centerOffset = (self.levels_n['left'][0][ind] + self.levels_n['right'][0][ind]) / 2.0

            x1, y1 = self.x_v[ind],     y_center[ind]
            x2, y2 = self.x_v[ind + 1], y_center[ind + 1]

            # Now just call afXShape with wbf
            y_temp, z_temp = self.afXShape(wbf, n=self.xshapePoints)

            y_temp += centerOffset

            real_x, real_y = functions.sn_to_xy(x1, y1, x2, y2, y_temp)
            if len(real_x) == 0:
                real_x = [x1] * self.xshapePoints
                real_y = [y1] * self.xshapePoints

            for i in range(min(len(real_x), self.xshapePoints)):
                xshape_lines[i].append((real_x[i], real_y[i], z_temp[i]))

            center_z.append(self.calCenter_z(real_x, real_y, z_temp, x1, y1))

        center_z.append(center_z[-1])

        # Flatten lines into out_x, out_y, out_z if needed
        for spoke in xshape_lines:
            for (xx, yy, zz) in spoke:
                out_x.append(xx)
                out_y.append(yy)
                out_z.append(zz)

        self.xshape_x = np.array(out_x)
        self.xshape_y = np.array(out_y)
        self.xshape_z = np.array(out_z)
        self.z_center = np.array(center_z)
        self.cross_section_type = 'AF'




    def getCenterlineElevation(self):
        '''Return z values for centerline.'''
        if self.xshape_x is None:
           self.setXShape()
        return self.z_center


#  ------to plot XShape---------
    def getXShapePlot(self):
        """
        Return matplotlib plot object showing the X-Shape of the channel.
        Now includes CF, PY, and AF cross-sections.
        """
        import numpy as np
        import matplotlib.pyplot as plt

        # Check if cross_section_type exists
        if hasattr(self, 'cross_section_type'):
            ctype = self.cross_section_type
        else:
            ctype = None  # Default to existing logic if undefined

        # Choose a middle station along the channel for plotting
        midInd = len(self.x_v) // 2

        # Compute total width at that station
        wbf = abs(self.levels_n["left"][0][midInd]) + abs(self.levels_n["right"][0][midInd])

        # 1) CF Cross-Section
        if ctype == 'CF':
            fig, ax = plt.subplots(1, 1)
            fig.suptitle('CF X-Shape for Channel')

            # Generate CF shape
            y, z = self.cfXShape(wbf)

            z = z + midInd * self.getPipeSlope()

            # Ensure the shape connects with banks
            y, z = self.addBankPoints(y, z, midInd)

            # Scale by dx if needed
            y = y * self.dx
            z = z * self.dx

            # Plot
            ax.plot(y, z, 'k-o', label=f'CF shape @ station {midInd}')
            ax.set_xlabel('Y (local cross-section)')
            ax.set_ylabel('Z')
            ax.legend()
            return fig

        # 2)PY Cross-Section
        if ctype == 'PY':
            fig, ax = plt.subplots(1, 1)
            fig.suptitle('PY X-Shape (Wavy U)')

            # Generate PY shape
            y, z = self.pyXShape(wbf)

            # Ensure the shape connects with banks
            y, z = self.addBankPoints(y, z, midInd)

            # Scale by dx if needed
            y = y * self.dx
            z = z * self.dx

            # Plot
            ax.plot(y, z, 'k-o', label=f'PY shape @ station {midInd}')
            ax.set_xlabel('Y')
            ax.set_ylabel('Z')
            ax.legend()
            return fig

        # 3) New AF (Angled Flat) Cross-Section
        if ctype == 'AF':
            fig, ax = plt.subplots(1, 1)
            fig.suptitle('AF X-Shape (Angled Flat)')

            # We call afXShape with the stored self.d1, etc. 
            y, z = self.afXShape(wbf, n=self.xshapePoints)

            # Ensure the shape connects with banks
            y, z = self.addBankPoints(y, z, midInd)

            # Possibly scale by dx
            y = y * self.dx
            z = z * self.dx

            ax.plot(y, z, 'b-o', label=f'AF shape @ station {midInd}')
            ax.set_xlabel('Y')
            ax.set_ylabel('Z')
            ax.legend()
            return fig

        # 4) AU/SU Cross-Sections
        cur_v = self.getDynamicCurv()
        maxCur = np.amax(cur_v)
        minCur = np.amin(cur_v)

        if maxCur == minCur or self.tz != -1:
            fig, ax = plt.subplots(1, 1)
            fig.suptitle('X-Shape for Channel')

            if self.tz == -1:
                y, z = self.pointXShape(midInd, maxCur, wbf, self.xshapePoints)
            else:
                y, z = self.suXShape(midInd, wbf, self.tz, self.xshapePoints)

            z = z + midInd * self.x_slope
            y, z = self.addBankPoints(y, z, midInd)

            y = y * self.dx
            z = z * self.dx

            ax.plot(y, z, 'k-', marker='o', label='x = ' + str(midInd))
            ax.set_xlabel('Y (related to center of channel)')
            ax.set_ylabel('Z')
            ax.legend()
            return fig

        # 5) Max Curvature vs. Min Curvature
        abs_cur_v = np.absolute(cur_v)
        fig, ax = plt.subplots(2, 1, sharex=True)
        fig.suptitle('Max Curvature X-Shape vs. Zero Curvature X-Shape')

        # Subplot #1: Max Curvature
        plt.subplot(212)
        indMax = np.argmax(abs_cur_v)
        wbf = abs(self.levels_n["left"][0][indMax]) + abs(self.levels_n["right"][0][indMax])
        y, z = self.pointXShape(indMax, maxCur, wbf, self.xshapePoints)
        si = self.getCenterline_sn()[indMax]
        z = z + si * self.getPipeSlope()
        y, z = self.addBankPoints(y, z, indMax)
        plt.plot(y, z, 'k-', marker='o', label=f'Max Curvature: x={indMax}')
        plt.xlabel('Y (related to center of channel)')
        plt.ylabel('Z')
        plt.legend()

        # Subplot #2: Min Curvature
        plt.subplot(211)
        indMin = np.argmin(abs_cur_v)
        wbf = abs(self.levels_n["left"][0][indMin]) + abs(self.levels_n["right"][0][indMin])
        y, z = self.pointXShape(indMin, minCur, wbf, self.xshapePoints)
        si = self.getCenterline_sn()[indMin]
        z = z + si * self.getPipeSlope()
        y, z = self.addBankPoints(y, z, indMin)

        y = y * self.dx
        z = z * self.dx

        plt.plot(y, z, 'k-', marker='o', label=f'Min Curvature: x={indMin}')
        plt.ylabel('Z')
        plt.legend()
        return fig



    def setXShape(self, n=-1):
        '''Calculate x, y, z values for Xshape of the whole channel.
           Also calculate the z values of centerline.
           xshapePointsDict: {(x, y): [z, (x_center, y_center)]}
        '''
        out_x, out_y, out_z = [], [], []
        xshapePointsList = []
        center_z = []
        y_center = self.getCenterline_y()
        s_center = self.getCenterline_sn()
        pipe_slope = self.getPipeSlope()

        cur_v = self.getDynamicCurv()
        maxCur = np.amax(np.absolute(cur_v))

        asFlag = True  # asymmetric flag
        if n != -1:
            asFlag = False  # innerPipePoints dict will be empty

        xshape_lines = [[] for i in range(self.xshapePoints)]

        for ind in range(len(y_center)-1):
            wbf = abs(self.levels_n['left'][0][ind]) + abs(self.levels_n['right'][0][ind])
            centerOffset = (self.levels_n['left'][0][ind] + self.levels_n['right'][0][ind])/2
            x1 = self.x_v[ind]
            y1 = y_center[ind]

            x2 = self.x_v[ind+1]
            y2 = y_center[ind+1]

            s = s_center[ind]
            # This if statement will determine whether it is AU or SU
            if asFlag:
                y_temp, z = self.pointXShape(ind, maxCur, wbf, self.xshapePoints)
            else:
                y_temp, z = self.suXShape(ind, wbf, n, self.xshapePoints)

            y_temp = y_temp + centerOffset
            real_x, real_y = functions.sn_to_xy(x1, y1, x2, y2, y_temp)
            # the following line may need to be commented out
            #z = z - pipe_slope*s                                            # use s instead of x
###############################################
#            if asFlag:
            for i in range(len(xshape_lines)):
                xshape_lines[i].append((real_x[i], real_y[i], z[i]))
#            else:
#                out_x += real_x.tolist()
#                out_y += real_y.tolist()
#                out_z += z.tolist()

            #find z for center line
            center_z.append(self.calCenter_z(real_x, real_y, z, x1, y1))

        center_z.append(center_z[-1])

        x_min = floor(min(self.innerPipePoints.keys()))
        x_max = ceil(max(self.innerPipePoints.keys())) + 1
        markPoints = [[] for i in range(ceil(x_max) - min(floor(x_min), 0))]

        for i in range(len(self.levels_x['left'][0])):
            x = self.levels_x['left'][0][i]
            y = self.levels_y['left'][0][i]
            z = self.levels_z['left'][0][i]
            markPoints[int(x)].append((y, z))

        for i in range(len(self.levels_x['right'][0])):
            x = self.levels_x['right'][0][i]
            y = self.levels_y['right'][0][i]
            z = self.levels_z['right'][0][i]
            markPoints[int(x)].append((y, z))

        for line in xshape_lines:
            line = functions.deleteCycles(line)
            for (x, y, z) in line:
                markPoints[x].append((y, z))

        for x in self.innerPipePoints.keys():
            innerPoint_y = self.innerPipePoints[x]
            xshape_yz = markPoints[x]
            if len(xshape_yz) == 0:
                continue
            xshape_yz.sort()
            xshape_y = [y for (y, z) in xshape_yz]
            xshape_z = [z for (y, z) in xshape_yz]
            for y in innerPoint_y:
                ind1, ind2 = functions.indexBound(y, xshape_y)
                if ind1 == ind2:
                    z = xshape_z[ind1]
                else:
                    z1 = xshape_z[ind1]
                    z2 = xshape_z[ind2]
                    y1 = xshape_y[ind1]
                    y2 = xshape_y[ind2]
                    alpha = (y-y1)/(y2-y1)
                    z = z1*(1-alpha) + z2*alpha
                out_x.append(x)
                out_y.append(y)
                out_z.append(z)

        self.xshape_x = np.array(out_x)
        self.xshape_y = np.array(out_y)
        self.xshape_z = np.array(out_z)
        self.z_center = np.array(center_z)
    
    def setXShapeCF(self):
        """
        Calculate x, y, z for the custom CF shape: z = (a * |x - b|)^c.
        We then store them in self.xshape_x, self.xshape_y, self.xshape_z.
        Also calculate the centerline z's in self.z_center.
        """
        import numpy as np
        from math import floor, ceil

        a = 1024
        b = 0.5
        c = 2.0

        out_x, out_y, out_z = [], [], []
        center_z = []

        y_center = self.getCenterline_y()
        s_center = self.getCenterline_sn()
        pipe_slope = self.getPipeSlope()

        # We won't use dynamic curvature for CF; often best to set it to 0
        # or just ignore it:
        # cur_v = self.getDynamicCurv()
        # maxCur = np.amax(np.absolute(cur_v))

        # We'll keep the same number of cross-section "spokes" as we do in setXShape()
        xshape_lines = [[] for _ in range(self.xshapePoints)]

        for ind in range(len(y_center) - 1):
            # channel "width" (left+right)
            wbf = abs(self.levels_n['left'][0][ind]) + abs(self.levels_n['right'][0][ind])
            centerOffset = (self.levels_n['left'][0][ind] + self.levels_n['right'][0][ind]) / 2

            x1, y1 = self.x_v[ind],     y_center[ind]
            x2, y2 = self.x_v[ind + 1], y_center[ind + 1]

            s = s_center[ind]

            # 2.1) Get local y and z for this cross-section
            #     Here we call a small helper function (below) that returns
            #     y_temp (the local cross-stream coordinate) and z (the shape).
            y_temp, z = self.cfXShape(wbf)

            # shift the channel so that the center is at 'centerOffset'
            y_temp = y_temp + centerOffset

            # 2.2) Convert local coordinates to global XY
            real_x, real_y = functions.sn_to_xy(x1, y1, x2, y2, y_temp)

            # Optionally incorporate pipe slope if that is your standard approach
            # z = z - pipe_slope * s

            # 2.3) Collect lines: same pattern as setXShape
            for i in range(len(xshape_lines)):
                xshape_lines[i].append((real_x[i], real_y[i], z[i]))

            # 2.4) centerline z
            center_z.append(self.calCenter_z(real_x, real_y, z, x1, y1))

        # add the last point’s center_z
        center_z.append(center_z[-1])

        # 2.5) Combine with bank & inner-pipe points (same as setXShape)
        x_min = floor(min(self.innerPipePoints.keys()))
        x_max = ceil(max(self.innerPipePoints.keys())) + 1
        markPoints = [[] for _ in range(ceil(x_max) - min(floor(x_min), 0))]

        # left bank
        for i in range(len(self.levels_x['left'][0])):
            xx = self.levels_x['left'][0][i]
            yy = self.levels_y['left'][0][i]
            zz = self.levels_z['left'][0][i]
            markPoints[int(xx)].append((yy, zz))

        # right bank
        for i in range(len(self.levels_x['right'][0])):
            xx = self.levels_x['right'][0][i]
            yy = self.levels_y['right'][0][i]
            zz = self.levels_z['right'][0][i]
            markPoints[int(xx)].append((yy, zz))

        # cross-section lines
        for line in xshape_lines:
            line = functions.deleteCycles(line)
            for (xx, yy, zz) in line:
                markPoints[xx].append((yy, zz))

        # 2.6) Interpolate z for the "innerPipePoints"
        for xx in self.innerPipePoints.keys():
            innerPoint_y = self.innerPipePoints[xx]
            xshape_yz = markPoints[xx]
            if len(xshape_yz) == 0:
                continue
            xshape_yz.sort()
            xshape_y = [yy for (yy, zz) in xshape_yz]
            xshape_z = [zz for (yy, zz) in xshape_yz]
            for yval in innerPoint_y:
                ind1, ind2 = functions.indexBound(yval, xshape_y)
                if ind1 == ind2:
                    zz = xshape_z[ind1]
                else:
                    z1, z2 = xshape_z[ind1], xshape_z[ind2]
                    y1, y2 = xshape_y[ind1], xshape_y[ind2]
                    alpha = (yval - y1) / (y2 - y1)
                    zz = z1 * (1 - alpha) + z2 * alpha

                out_x.append(xx)
                out_y.append(yval)
                out_z.append(zz)

        self.xshape_x = np.array(out_x)
        self.xshape_y = np.array(out_y)
        self.xshape_z = np.array(out_z)
        self.z_center = np.array(center_z)

    def setXShapePY(self):
        """
        Build the geometry for a 'PY' cross-section, calling pyXShape(...) at each station
        and merging with bank points, just like your AU or CF methods.
        """


        out_x, out_y, out_z = [], [], []
        xshape_lines = [[] for _ in range(self.xshapePoints)]
        center_z = []

        y_center = self.getCenterline_y()
        s_center = self.getCenterline_sn()
        pipe_slope = self.getPipeSlope()

        # Loop over each station
        for ind in range(len(y_center)-1):
            wbf = abs(self.levels_n['left'][0][ind]) + abs(self.levels_n['right'][0][ind])
            centerOffset = (self.levels_n['left'][0][ind] + self.levels_n['right'][0][ind]) / 2

            x1, y1 = self.x_v[ind],     y_center[ind]
            x2, y2 = self.x_v[ind+1],   y_center[ind+1]

            s = s_center[ind]

            # 1) local cross-section for this station
            y_temp, z_temp = self.pyXShape(wbf)

            # shift cross-section horizontally by centerOffset
            y_temp = y_temp + centerOffset

            # 2) convert local y_temp => global coordinates
            real_x, real_y = functions.sn_to_xy(x1, y1, x2, y2, y_temp)

            # optionally slope the shape:
            # z_temp = z_temp - pipe_slope*s

            for i in range(len(xshape_lines)):
                xshape_lines[i].append((real_x[i], real_y[i], z_temp[i]))

            # centerline z
            center_z.append(self.calCenter_z(real_x, real_y, z_temp, x1, y1))

        # add last station center
        center_z.append(center_z[-1])

        # Merge cross-section lines with left/right bank points
        x_min = floor(min(self.innerPipePoints.keys()))
        x_max = ceil(max(self.innerPipePoints.keys())) + 1
        markPoints = [[] for _ in range(ceil(x_max) - min(floor(x_min), 0))]

        # banks: left + right
        for i in range(len(self.levels_x['left'][0])):
            xx = self.levels_x['left'][0][i]
            yy = self.levels_y['left'][0][i]
            zz = self.levels_z['left'][0][i]
            markPoints[int(xx)].append((yy, zz))

        for i in range(len(self.levels_x['right'][0])):
            xx = self.levels_x['right'][0][i]
            yy = self.levels_y['right'][0][i]
            zz = self.levels_z['right'][0][i]
            markPoints[int(xx)].append((yy, zz))

        # cross-section lines
        for line in xshape_lines:
            line = functions.deleteCycles(line)
            for (xx, yy, zz) in line:
                markPoints[xx].append((yy, zz))

        # fill in "innerPipePoints"
        for xx in self.innerPipePoints.keys():
            innerPoint_y = self.innerPipePoints[xx]
            xshape_yz = markPoints[xx]
            if len(xshape_yz) == 0:
                continue

            xshape_yz.sort()  # sort by Y
            xshape_y = [pair[0] for pair in xshape_yz]
            xshape_z = [pair[1] for pair in xshape_yz]

            for yval in innerPoint_y:
                ind1, ind2 = functions.indexBound(yval, xshape_y)
                if ind1 == ind2:
                    zz = xshape_z[ind1]
                else:
                    z1, z2 = xshape_z[ind1], xshape_z[ind2]
                    yy1, yy2 = xshape_y[ind1], xshape_y[ind2]
                    alpha = (yval - yy1) / (yy2 - yy1)
                    zz = z1*(1 - alpha) + z2*alpha

                out_x.append(xx)
                out_y.append(yval)
                out_z.append(zz)

        self.xshape_x = np.array(out_x)
        self.xshape_y = np.array(out_y)
        self.xshape_z = np.array(out_z)
        self.z_center = np.array(center_z)
        self.cross_section_type = 'PY'
    

    def addBoulders(self, num, size_mean, size_std, height):
        '''
        Add boulders 
        avail_pts - nested list;
                    elem: [set(available y values), (y1, z1), ...]
        '''
        size_mean = size_mean/self.dx
        size_std = size_std/self.dx
        height = height/self.dx

        x_min = np.amin(self.xshape_x)
        x_min = int(min(x_min, 0))
        x_max = int(np.amax(self.xshape_x) + 1)
        avail_pts = [[set()] for i in range(x_min, x_max)]

        for i in range(len(self.xshape_x)):
            x = int(self.xshape_x[i])
            avail_pts[x][0].add(self.xshape_y[i])

        area = []
        check_x = set(list(range(x_min, x_max)))
        
        while num > 0:
            area, check_x = self.cutArea(avail_pts, size_mean, size_std, check_x, x_min, x_max)
            if area == []:
                break
            boulder = self.createBoulder(area, height)
            self.updateBoulder(boulder)
            num -= 1


    def addCheckDam(self, loc, height, thick):
        '''
        Add check dam
        loc - location along meandering stream.
        height - height from the centerline point.
        thick - how thick is the dam.
        '''
        height = height/self.dx
        thick = thick/self.dx

        loc_ind = np.where(self.s_center > loc)[0]
        loc_ind = np.amin(loc_ind)
        s = self.getSlope()[loc_ind]
        x_cp = self.x_v[loc_ind]
        y_cp = self.y_center[loc_ind]
        z_cp = self.z_center[loc_ind]
        lf_range = np.amax(self.levels_n['left'][0])
        rt_range = np.amin(self.levels_n['right'][0])

#        x_len_inc, y_len_inc, x_wid_inc, y_wid_inc = 0, 0, 0, 0

        if abs(s) == math.inf:
            x_len_inc = 1
            y_len_inc = 0
            x_wid_inc = 0
            y_wid_inc = 1
        elif s == 0:
            x_len_inc = 0
            y_len_inc = 1
            x_wid_inc = 1
            y_wid_inc = 0
        elif abs(s) > 1:
            x_len_inc = 1
            y_len_inc = -1/s
            x_wid_inc = 1/s
            y_wid_inc = 1
        else:
            x_len_inc = s
            y_len_inc = -1
            x_wid_inc = 1
            y_wid_inc = 1/s

        pt_crt_x, pt_crt_y = round(x_cp), round(y_cp)
        ck_dam_pts = []

        for dummy in range(int(lf_range)):
            ck_dam_pts.append((pt_crt_x, pt_crt_y))
            pt_crt_x = round(pt_crt_x - x_len_inc)
            pt_crt_y = round(pt_crt_y - y_len_inc)

            for i in range(thick):
                pt_wid_x = round(pt_crt_x + i*x_wid_inc)
                pt_wid_y = round(pt_crt_y + i*y_wid_inc)
                ck_dam_pts.append((pt_wid_x, pt_wid_y))

        pt_crt_x, pt_crt_y = round(x_cp), round(y_cp)
        for dummy in range(abs(int(rt_range))):
            ck_dam_pts.append((pt_crt_x, pt_crt_y))
            pt_crt_x = round(pt_crt_x + x_len_inc)
            pt_crt_y = round(pt_crt_y + y_len_inc)

            for i in range(thick):
                pt_wid_x = round(pt_crt_x + i*x_wid_inc)
                pt_wid_y = round(pt_crt_y + i*y_wid_inc)
                ck_dam_pts.append((pt_wid_x, pt_wid_y))

        for (x, y) in ck_dam_pts:
            ind_x = np.where(self.xshape_x == x)[0]
            ind_y = np.where(self.xshape_y == y)[0] 
#                if len(np.intersect1d(ind_x, ind_y)) == 0:
#                    print('ind_x', ind_x)
#                    print('ind_y', ind_y)
#                    print('x, y', x, y)
            inter = np.intersect1d(ind_x, ind_y)
            if len(inter) > 0:
                ind = np.intersect1d(ind_x, ind_y)[0]
                self.xshape_z[ind] = z_cp + height

    def tolist(self):
        '''Return x, y, z values for all levels in a secondary list'''
        x = []
        y = []
        z = []

        x += self.x_v.tolist()
        y += self.getCenterline_y().tolist()
        z += self.getThalweg().tolist()

        self.helpAppend(x, self.levels_x)
        self.helpAppend(y, self.levels_y)
        self.helpAppend(z, self.levels_z)

        return [x, y, z]


    def tocsv(self, outfile):
        '''Outwrite xy values and xz values of all banks to output file.'''
        header = ["x", "y", 'z']
        out = [header]
        xyz = self.tolist()
        xyz_out = [[round(xyz[0][i]*self.dx, 3), round(xyz[1][i]*self.dx, 3), round(xyz[2][i]*self.dx, 3)] for i in range(len(xyz[0]))]
        out += xyz_out
        with open(outfile+".csv", 'w') as cf:
            cw = csv.writer(cf)
            cw.writerows(out)


    def __str__(self):
        '''Turn most important data of river.'''
        sl = self.getSL()
        aveWbf = self.getAveWbf()
        aveHbf = self.getAveHbf()
        slope = self.getRiverSlope()
        coWbf = self.getCoWbf()
        coHbf = self.getCoHbf()
        s = 'Sinuosity:'+str(round(sl, 3))+'\n'
        s += 'Channel Slope:'+str(slope)+'\n'
        s += 'Average Width of Inner Channel:'+str(round(aveWbf, 3))+'\n'
        s += 'Average Height of Inner Channel:'+str(round(aveHbf, 3))+'\n'
        s += 'Coefficient of Variation (W_ic):'+str(round(coWbf, 3))+'\n'
        s += 'Coefficient of Variation (H_ic):'+str(round(coHbf, 3))+'\n'

        if len(self.levels_n['left']) == 1:
            return s
        for i in range(1, len(self.levels_n['left'])):
            s += 'Average Width Offset of '+'L'+str(i)+' Outer Bank is: '+str(round(np.average(self.levels_n['left'][i])*self.dx, 3)) + '\n'
            
        if len(self.levels_n['right']) == 1:
            return s
        for i in range(1, len(self.levels_n['right'])):
            s += 'Average Width Offset of '+'R'+str(i)+' Outer Bank is: '+str(abs(round(np.average(self.levels_n['right'][i])*self.dx, 3))) + '\n'

        return s


    def constructChannel(self):
        ''' Construct channel based on the information stored in Channel.

        The construction follows following steps:
            1. Build up centerline.
            2. Build up thalweg.
            3. Build up inner channels.
            4. Build up xshape points.
            5. Build up outer banks.
        '''
        # Build up centerline.
        self.setCenterline()
        self.setCenterline_sn()
        self.loopCenterline(goal)

        # Build up thalweg.
        self.setThalweg()

        # Build up inner channels.
        self.createInnerChannel()
        self.loopHbf(goal)


##############################################################

# hELPER FUNCTIONS

    def helpAppend(self, li, dic):
        for array in dic["left"]:
            li += array.tolist()
        for array in dic["right"]:
            li += array.tolist()


    def pointXShape(self, ind, maxCur, wbf, n):
        '''Return y values and z values of XSection of given x

        n -- number of points to calculate XSection
        '''
        cur = self.getDynamicCurv()[ind]
        pipe_slope = self.getPipeSlope()
        si = self.getCenterline_sn()[ind]
        xVal = np.round(self.x_v[ind])

        if maxCur == 0:
            B = 1/2
        else:
            B = 1/2 * (1 - abs(cur/maxCur))

        if B < 0.1:
            B = 0.1

        if B == 1:
            L = 1
        else:
            L = -1*log(2)/log(B)

        lbx = self.levels_x['left'][0]
        lb = self.levels_z['left'][0]
        rbx = self.levels_x['right'][0]
        rb = self.levels_z['right'][0]
        lb_ind = np.where(lbx == xVal)[0]
        if len(lb_ind) == 0:
            bankH = rb[np.where(rbx == xVal)[0][0]]
        else:
            bankH = lb[lb_ind[0]]
        hbf = bankH - self.thalweg[ind]

        n_y = np.array([-wbf*x/(n+1) + (1/2)*wbf  for x in range(1, n+1)])
        Y = (wbf/2-n_y) / wbf

        if cur <= 0:
            n_z = 4 * hbf * (Y**L) * (1-Y**L)
        else:
            n_z = 4 * hbf * ((1-Y)**L) * (1-(1-Y)**L)

        n_z = self.thalweg[ind] + hbf - n_z
        return n_y, n_z


    def suXShape(self, ind, wbf, tzn, n):
        '''Return y values and z values of Symmetric XSection of given x.
        
        ind - index of x value on centerline
        tzn - number of points on base
        n - number of points on XS
        '''
        xVal = np.round(self.x_v[ind])
        lbx = self.levels_x['left'][0]
        lb = self.levels_z['left'][0]
        rbx = self.levels_x['right'][0]
        rb = self.levels_z['right'][0]
        lb_ind = np.where(lbx == xVal)[0]
        if len(lb_ind) == 0:
            bankH = rb[np.where(rbx == xVal)[0][0]]
        else:
            bankH = lb[lb_ind[0]]
        hbf = bankH - self.thalweg[ind]
        n_y = np.array([-wbf*x/(n+1) + (1/2)*wbf  for x in range(1, n+1)])

        n_z = []
        sidePoints = floor(((n-tzn)/2) + 1)

        if (n-tzn) % 2 != 0:
            tzn += 1

        for i in range(1, sidePoints):
            n_z.append((i/sidePoints)*hbf)

        zEnd = n_z[::-1]
        n_z += [hbf] * tzn
        n_z += zEnd

        n_z = self.thalweg[ind] + hbf - n_z
        return n_y, np.array(n_z)

    def cfXShape(self, wbf, hbf=None, thalweg=None):
        """
        Returns (y, z) for a symmetrical cross‐section using a scaled parametric
        equation. The profile is designed so that:
          - at the center (y=0): elevation = thalweg (lowest point)
          - at the banks (y = ±wbf/2): elevation = thalweg + hbf (bank level)
          
        If hbf or thalweg are not provided, they are computed using the mid–station.
        """
        if hbf is None or thalweg is None:
            midInd = len(self.x_v) // 2
            xVal = np.round(self.x_v[midInd])
            lbx = self.levels_x['left'][0]
            lb = self.levels_z['left'][0]
            rbx = self.levels_x['right'][0]
            rb = self.levels_z['right'][0]
            lb_ind = np.where(lbx == xVal)[0]
            if len(lb_ind) == 0:
                bankH = rb[np.where(rbx == xVal)[0][0]]
            else:
                bankH = lb[lb_ind[0]]
            hbf = bankH - self.thalweg[midInd]
            thalweg = self.thalweg[midInd]
        
        # Create local y coordinates spanning the full width:
        y_temp = np.linspace(-wbf/2, wbf/2, self.xshapePoints)
        # A normalized parabolic (U–shaped) profile:
        a = 1.0
        b = 2.0
        c = 5.0
        profile = a - (b* np.abs(y_temp) / wbf)**c
        # profile = a*(np.abs(y_temp - b)**c)  
        # Map profile to elevation:
        #   When profile=1, elevation = thalweg;
        #   When profile=0, elevation = thalweg + hbf.
        z_temp = thalweg + hbf * (1 - profile)  
    
        return y_temp, z_temp

    def pyXShape(self, wbf, hbf=None, thalweg=None, A=10.0, B=3.0, freq=2):
        """
        Returns (y, z) arrays for a 'wavy U–shape' cross–section based on:
            z = A*(1 - x^2) + B*sin(2*pi*freq*x),
        where x is a dimensionless coordinate spanning [-1, 1] and y is mapped to
        [-wbf/2, wbf/2]. The raw z profile is then scaled so that its maximum maps to
        thalweg (lowest elevation) and its minimum maps to thalweg+hbf (bank level),
        effectively flipping the profile upside down relative to the original.
        
        If hbf or thalweg are not provided, they are computed using the mid–station.
        User-provided values for A, B, and freq control the shape.
        """
        # If hbf or thalweg are not provided, compute them at the mid–station:
        if hbf is None or thalweg is None:
            midInd = len(self.x_v) // 2
            xVal = np.round(self.x_v[midInd])
            lbx = self.levels_x['left'][0]
            lb = self.levels_z['left'][0]
            rbx = self.levels_x['right'][0]
            rb = self.levels_z['right'][0]
            lb_ind = np.where(lbx == xVal)[0]
            if len(lb_ind) == 0:
                bankH = rb[np.where(rbx == xVal)[0][0]]
            else:
                bankH = lb[lb_ind[0]]
            hbf = bankH - self.thalweg[midInd]
            thalweg = self.thalweg[midInd]
        
        n = self.xshapePoints  # number of sample points across the channel
        # Dimensionless cross-stream coordinate:
        x_vals = np.linspace(-1, 1, n)
        # Map x to local y coordinates spanning [-wbf/2, wbf/2]:
        y_vals = x_vals * (wbf / 2.0)
    
        # Compute the raw profile using the given parameters:
        z_raw = A * (1 - x_vals**2) + B * np.sin(2 * np.pi * freq * x_vals)
    
        # Determine the range of the raw profile:
        z_min = np.min(z_raw)
        z_max = np.max(z_raw)
    
        # Reverse the mapping: when z_raw is at its minimum, we want elevation = thalweg+hbf;
        # when z_raw is at its maximum, we want elevation = thalweg.
        z_vals = thalweg + hbf * (1 - (z_raw - z_min) / (z_max - z_min))
    
        return y_vals, z_vals

    def afXShape(self, wbf, n=21, hbf=None, thalweg=None):
        """
        Build the trapezoid using only the stored self.d1, self.d2, self.ang1, self.ang2.
        """
        import numpy as np, math

        # If hbf or thalweg not provided, get from mid-station
        midInd = len(self.x_v) // 2
        if hbf is None:
            hbf = self.getHbf()
        if thalweg is None:
            thalweg = self.getThalweg()[midInd]

        # Use the values stored in self
        d1   = self.d1
        d2   = self.d2
        ang1 = math.radians(self.ang1)
        ang2 = math.radians(self.ang2)

        w1 = d1 / math.tan(ang1) if math.tan(ang1) != 0 else 0
        w2 = d2 / math.tan(ang2) if math.tan(ang2) != 0 else 0

        remainingWidth = wbf - (w1 + w2)
        if remainingWidth < 0:
            remainingWidth = 0

        leftBank  = -wbf / 2.0
        leftBed   = leftBank + w1
        rightBed  = leftBed + remainingWidth
        rightBank = rightBed + w2

        total_length = wbf
        np_left = max(2, int(round(n * (w1 / total_length))))
        np_flat = max(2, int(round(n * (remainingWidth / total_length)))) if remainingWidth > 0 else 0
        np_right = max(2, int(round(n * (w2 / total_length))))
        total_pts = np_left + np_flat + np_right
        if total_pts < n:
            np_flat += (n - total_pts)

        left_x = np.linspace(leftBank, leftBed, num=np_left, endpoint=True)
        left_z = np.linspace(0, -d1, num=np_left, endpoint=True)

        if np_flat > 0:
            flat_x = np.linspace(leftBed, rightBed, num=np_flat, endpoint=True)
            flat_z = np.full_like(flat_x, -max(d1, d2))
        else:
            flat_x = np.array([])
            flat_z = np.array([])

        right_x = np.linspace(rightBed, rightBank, num=np_right, endpoint=True)
        right_z = np.linspace(-d2, 0, num=np_right, endpoint=True)

        y_vals = np.concatenate([left_x, flat_x, right_x])
        z_raw  = np.concatenate([left_z, flat_z, right_z])

        top = thalweg + hbf
        bed_elev = top - max(d1, d2)
        norm = (0 - z_raw) / max(d1, d2) if max(d1, d2) != 0 else 0
        z_vals = top - norm * (top - bed_elev)

        return y_vals, z_vals




    def addBankPoints(self, y, z, ind):
        '''Add bank points to xshape points.
        
        y - y values for xshape points
        z - z values for xshape points
        ind - where the xshape points are calculated

        Return:
        y, z with bank points added
        '''
        leftEdge = y[0]-(y[1]-y[0])
        y = np.append(leftEdge, y)
        rightEdge = y[-1]+(y[1]-y[0])
        y = np.append(y, rightEdge)

        z = np.append(self.levels_z['left'][0][0], z)
        z = np.append(z, self.levels_z['left'][0][0])

        for i in range(1, len(self.levels_n['left'])):
            y = np.append(self.levels_n['left'][i][ind] - self.levels_n['left'][0][ind] + leftEdge, y)
            z = np.append(self.levels_z['left'][i][0], z)


        for i in range(1, len(self.levels_n['right'])):
            y = np.append(y, self.levels_n['right'][i][ind] - self.levels_n['right'][0][ind] + rightEdge)
            z = np.append(z, self.levels_z['right'][i][0])
        return y, z


    def setDynamicCurv(self):
        ''' Calculate the dynamic curve of the centerline.
        '''
        x_v = self.getCenterline_x()
        y_v = self.getCenterline_y()
        slopeVectorList = []

        for i in range(len(x_v)-1):
            v = (x_v[i+1]-x_v[i], y_v[i+1]-y_v[i])
            slopeVectorList.append(v)

        cur = []
        piCheck = pi/2
        for i in range(len(slopeVectorList)-1):
            v1 = slopeVectorList[i]
            v2 = slopeVectorList[i+1]
            angle = functions.angle_between(v1, v2)
            if np.cross(v1, v2) >= 0:
                cur.append(functions.angle_between(v1, v2))
            else:
                cur.append(functions.angle_between(v1, v2) * -1)

        cur.append(cur[-1])
        cur.append(cur[-1])
        self.dynamicCurv = np.array(cur)


    def calCenter_z(self, real_x, real_y, z, x1, y1):
        '''Calculate the z value for the centerline point.'''
        xpoints = [(real_x[i], real_y[i]) for i in range(len(real_x))]
        minInd = 0
        minDist = functions.pointDist(xpoints[0], (x1, y1))
        diff = [minDist]
        for i in range(1, len(xpoints)):
            dist = functions.pointDist(xpoints[i], (x1, y1))
            diff.append(dist)
            if dist < minDist:
                minInd = i
                minDist = dist

        if minDist == 0 or (minInd-1 < 0 and minInd+1 >= len(diff)) :
            return z[minInd]
        elif minInd+1 < len(diff) and (minInd-1 < 0 or diff[minInd-1] >= diff[minInd+1]):
            minInd2 = minInd+1
        else:
            minInd2 = minInd-1
        
        z1 = z[minInd]
        z2 = z[minInd2]
        minX1 = real_x[minInd]
        minX2 = real_x[minInd2]
        if minX1 == minX2:
            return (z1+z2)/2
        elif minX1 < minX2:
            alpha = (x1-minX1)/(minX2-minX1)
            return alpha*z1 + (1-alpha)*z2
        else:
            alpha = (x1-minX2)/(minX1-minX2)
            return alpha*z1 + (1-alpha)*z2


    def updateXShapePointsList(self, pointsList, x_v, y_v, z_v, x_center, y_center):
        '''Update the XShape points in XShape points Dictionary.

        pointsList -- [(x, y, z)]
        '''
##################################
#        checkRange = max(int(np.amax(y_v) - np.amin(y_v)), int(np.amax(x_v) - np.amin(x_v)))
#        x_check = self.getCenterline_x()
#        y_check = self.getCenterline_y()
#        inCount = outCount = 0
#        for i in range(len(x_v)):
#            x = int(x_v[i])
#            y = int(y_v[i])
#
#            left = max(int(x_center - checkRange), 0)
#            right = min(int(x_center + checkRange), len(x_check))
#
#            dist = np.square(x_check[left:right] - x) + np.square(y_check[left:right] - y)
#            minIndex = np.argmin(dist)
#            minDistX = int(x_check[left:right][minIndex])
#
#            if minDistX == int(x_center):
#                pointsList.append((x, y, z_v[i]))
##################################
        lbx = self.levels_x['left'][0]
        lby = self.levels_y['left'][0]
        rbx = self.levels_x['right'][0]
        rby = self.levels_y['right'][0]

        head = (x_v[0], y_v[0])
        tail = (x_v[-1], y_v[-1])

        head_dist_lb = np.square(lbx - head[0]) + np.square(lby - head[1])
        head_dist_rb = np.square(rbx - head[0]) + np.square(rby - head[1])

        tail_dist_lb = np.square(lbx - tail[0]) + np.square(lby - tail[1])
        tail_dist_rb = np.square(rbx - tail[0]) + np.square(rby - tail[1])

        if min(np.amin(head_dist_lb), np.amin(head_dist_rb)) < 10:
            for i in range(round(len(x_v)/2)):
                pointsList.append((x_v[i], y_v[i], z_v[i]))
        else:
            for i in range(round(len(x_v)/2)-4, round(len(x_v)/2)):
                pointsList.append((x_v[i], y_v[i], z_v[i]))


        if min(np.amin(tail_dist_lb), np.amin(tail_dist_rb)) < 10:
            for i in range(round(len(x_v)/2), len(x_v)):
                pointsList.append((x_v[i], y_v[i], z_v[i]))
        else:
            for i in range(round(len(x_v)/2), round(len(x_v)/2)+4):
                pointsList.append((x_v[i], y_v[i], z_v[i]))

##################################

        return pointsList


    def cutArea(self, avail_pts, size_mean, size_std, check, x_min, x_max):
        '''
        avail_pts - nested list;
                    elem: [set(available y values), (y1, z1), ...]
        '''
        find = False
        area = []
        length = int(np.random.normal(size_mean, size_std))
        length = max(length, 5)
        length = min(length, size_mean+3*size_std)
        width = int(np.random.normal(size_mean, size_std))
        width = max(width, 5)
        width = min(width, size_mean+3*size_std)
        width_half = width/2

        while not find:
            # check if no space left
            if len(check) == 0:
                break

            # pop out a random x position
            ind = random.sample(check, 1)[0]
            check.remove(ind)
            if ind - round(width_half) < x_min or \
                    ind + round(width_half) >= x_max:
                continue
            
            # all y values available at this x position
            y_pool = avail_pts[ind][0].copy()
            if len(y_pool) < length:
                continue
            
            # find a y that is valid
            check_y_pool = 0       # This check if y pool it self has valid ys.
            all_ok = True
            y_start = 0
            while len(y_pool) != 0:
                # pop out a random starting y value
                y_start = random.sample(y_pool, 1)[0]
                y_pool.remove(y_start)
                ys = list(range(y_start, y_start+length))
                all_ok = True

                # check x by x if there are enough space 
                for i in range(ind-round(width_half), ind+round(width_half)+1):
                    y_set = avail_pts[i][0]

                    for y in ys:
                        if y not in y_set:
                            all_ok = False
                            if i == ind:
                                check_y_pool += 1
                            break

#                    if not ok:
#                        break

                if all_ok:
                    break

            if check_y_pool < length:
                check.add(ind)

            if all_ok:
                for i in range(ind-round(width_half), ind+round(width_half)+1):
#                    area.append([])
                    for t in range(length):
                        avail_pts[i][0].remove(t+y_start)
#                        area[-1].append((i, t+y_start))
                area = [(ind-round(width_half), ind+round(width_half)), (y_start, y_start+length-1)]
                find = True
            
        return area, check


    def createBoulder(self, area, height=5):
        '''
        area - list of range
                [(x_start, x_end), (y_start, y_end)]
        '''
        # get length
        end_y, start_y = area[1][1], area[1][0]
        end_x, start_x = area[0][1], area[0][0]

        length = max(end_y-start_y+1, end_x-start_x+1)
        r = length/2
        temp_x = np.arange(length)
        temp_y = np.arange(length)
        temp_x, temp_y = np.meshgrid(temp_x, temp_y)
        z = np.sqrt(np.square(temp_x-r) + np.square(temp_y-r))
        z = (r-z)/r
        z = np.maximum(z, 0)
        
        for i in range(len(z)):
            for t in range(len(z[0])):
                if z[i][t] > 0.5:
                    z[i][t] = 0.5 + (z[i][t]-0.5)/2

        z = z/0.75 * height
        err = np.random.random_sample(z.shape)*height/10
        z = z + err
        
        diff_x = length - (end_x - start_x + 1)
        diff_y = length - (end_y - start_y + 1)
        if diff_x > 0:
            z = z[:, floor(diff_x/2):len(z[i])-ceil(diff_x/2)]
            temp_x = temp_x[:, floor(diff_x/2):len(temp_x[i])-ceil(diff_x/2)]
            temp_y = temp_y[:, floor(diff_x/2):len(temp_y[i])-ceil(diff_x/2)]
        elif diff_y > 0:
            z = z[floor(diff_y/2):len(z)-ceil(diff_y/2), :]
            temp_x = temp_x[floor(diff_y/2):len(temp_y)-ceil(diff_y/2), :]
            temp_y = temp_y[floor(diff_y/2):len(temp_y)-ceil(diff_y/2), :]

        start_x -= temp_x[0][0]
        start_y -= temp_y[0][0]
        x = temp_x + start_x
        y = temp_y + start_y

        return (x, y, z)


    def updateBoulder(self, boulder):
        (b_x, b_y, b_z) = boulder
        for i in range(len(b_x)):
            for t in range(len(b_x[0])):
                x = b_x[i][t]
                y = b_y[i][t]
                z = b_z[i][t]
                
                ind_x = np.where(self.xshape_x == x)[0]
                ind_y = np.where(self.xshape_y == y)[0] 
#                if len(np.intersect1d(ind_x, ind_y)) == 0:
#                    print('ind_x', ind_x)
#                    print('ind_y', ind_y)
#                    print('x, y', x, y)
                ind = np.intersect1d(ind_x, ind_y)[0]
                self.xshape_z[ind] += z


    def perlinThalweg(self, height):
        '''
        2D perlin function through the whole inner channel.
        height is the maximum difference of the noise.
        '''
        # decide the x range and y range of the channel base
        height = height/self.dx

        min_x = np.amin(self.xshape_x)
        min_y = np.amin(self.xshape_y)
        max_x = np.amax(self.xshape_x)
        max_y = np.amax(self.xshape_y)

        diff_x = max_x - min_x + 1
        diff_y = max_y - min_y + 1

        # create a frame for perlin2D function
        diff_x = ceil(diff_x/10)
        diff_y = ceil(diff_y/10)

        lin_x = np.linspace(0, diff_x, diff_x*10)
        lin_y = np.linspace(0, diff_y, diff_y*10)
        x, y = np.meshgrid(lin_x, lin_y)

        # generage 2d perlin noise
        z = functions.perlin2D(x, y)
        z *= height

        # update noise to channel base
        for i in range(len(self.xshape_x)):
            xi = self.xshape_x[i]
            yi = self.xshape_y[i]
            xi = xi - min_x
            yi = yi - min_y
            zi = z[yi][xi]

            self.xshape_z[i] += zi

