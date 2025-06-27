'''This is a basic pipe module.

A pipe consists of a centerline, unlimited number of levels 
on each side. The levels are first calculated in s-n coordinate 
system, and then they will be converted into x-y coordinate 
system.
'''

import numpy as np
from . import functions
from math import pi, sqrt, log, floor, ceil
import csv

class Pipe:

    def __init__(self, x_len=100, x_slope=0.01, dx=1, zd=1000):
        '''Pipe class initiator

        x_len -- int; length in x direction
        x_slope -- float; slope in x direction
        dx -- int; resolution in x direction, default is 1.
                if dx=0.1, it means the size of grid in innerPipe is 0.1.
        class private variables:
        x_v -- array; x values of centerline
        x_index -- array; index of x_v
        y_center -- array; y values of centerline
        levels_x -- dictionary;
                    key: string; directrion ('left' or 'right')
                    value: list of arrays of x values of levels
        levels_y -- dictionary;
                    key: string; directrion ('left' or 'right')
                    value: list of arrays of y values of levels
        levels_z -- dictionary;
                    key: string; directrion ('left' or 'right')
                    value: list of arrays of z values of levels
        levels_n -- dictionary;
                    key: string; directrion ('left' or 'right')
                    value: list of arrays of n values of levels
        leveldict -- dictionary;
                    key: direction of levels;
                    value: dict of levels of given direction
        s_center -- array; s values of centerline (sn coordinate)
        slope_center -- array; slope of centerline
        innerPipePoints -- dictionary;
                            key - x value; value - y values set
        '''
        self.x_len = x_len
        self.dx = dx
        self.x_slope = x_slope
        self.zd = zd/self.dx
        self.x_v = None
        self.x_index = np.array(list(range(ceil(x_len/dx))))
        self.y_center = None
        self.levels_x={"left":[], "right":[]}
        self.levels_y={"left":[], "right":[]}
        self.levels_z={"left":[], "right":[]}
        self.levels_n={"left":[], "right":[]}
        self.bounds = {'left': None, 'right': None}
        self.leveldict = {'x':self.levels_x,
                          'y':self.levels_y,
                          'z':self.levels_z,
                          'n':self.levels_n
                          }
        self.thalweg = None
        self.s_center = None
        self.slope_center = None
        self.innerPipePoints = {}


    def setCenterline(self, fun=None):
        '''Calculate y values of centerline

        fun -- function used to calculate centerline
        '''
        if fun is None:
            def fun(x):
                return x, np.array([0]*len(x))
        self.x_v, self.y_center = self.createCenterline(fun)


    def setCenterlineManual(self, y_v, x_v=None):
        ''' Manually set the centerline to y_v.
            If x_v is None, the original x_len will be used to
                calculate the x_v.
            If x_v is not None, it will replace the original x_v.
            The length of y_v will be cut or extend to match the length of x_v.

            y_v - numpy array; y values of the centerline
            x_v - numpy array; x values of the centerline
        '''
        if x_v is not None:
            self.x_v = x_v
        else:
            def fun(x):
                return x, np.array([0]*len(x))
            self.x_v, dummy = self.createCenterline(fun)

        diff = len(y_v) - len(self.x_v)
        if diff >= 0:
            self.y_center = y_v[:len(self.x_v)]
        else:
            self.y_center = np.concatenate((y_v, np.array([y_v[-1]]*(-diff))))


    def smoothCenterline(self, degree=1):
        ''' Smoothen centerline by taking average of neighbors.

        Input:
        degree - number of neighbors to take average 
                    (1 means 1 left neighbor and 1 rightneighbor,
                     so 3 points in total.)
        '''
        if degree == 0:
            return

        y_v = self.getCenterline_y()
        out = np.array([0.0]*len(y_v))
        for i in range(len(y_v)):
            left = i-degree
            right = i+degree+1
            if left < 0:
                left = 0
            if right > len(y_v):
                right = len(y_v)
            out[i] = np.average(y_v[left:right])
        self.y_center = out


    def getCenterline_y(self):
        '''Return self.y_center'''
        if self.y_center is None:
            self.setCenterline()
        return self.y_center


    def getCenterline_x(self):
        '''Return self.x_v'''
        if self.x_v is None:
            self.setCenterline()
        return self.x_v


    def setSlope(self):
        '''Calsulate slope for centerline.'''
        y_center = self.getCenterline_y()
        self.slope_center = functions.slope_v(self.x_v, y_center)


    def getSlope(self):
        '''Return self.slope_center.'''
        if self.slope_center is None:
            self.setSlope()
        return self.slope_center


    def setCenterline_sn(self):
        '''Calculate s values for centerline.'''
        y_v = self.getCenterline_y()
        self.s_center = self.calculate_s(self.x_v, y_v)


    def getCenterline_sn(self):
        '''Return self.s_center
        '''
        if self.s_center is None:
            self.setCenterline_sn()
        return self.s_center


    def setThalweg(self, fun=None, zd=0):
        '''Calculate z values for thalweg

        zd -- Datum
        fun -- function used to calculate thalweg

        Value will be modified:
        self.thalweg
        '''
        s_v = self.getCenterline_sn()

        if fun == None:
            x = s_v
            y = np.zeros(len(x))
        else:
            x, y = fun(s_v)
            y = y / self.dx
        dynamicSlope = self.getDynamicPipeSlope()
        delta_s_v = np.concatenate((np.array([0]), s_v[1:] - s_v[0:-1]))
        delta_z = np.multiply(dynamicSlope, delta_s_v)
        z_v = np.empty_like(delta_z)
        z_v[0] = 0
        for i in range(1, len(z_v)):
            z_v[i] = z_v[i-1] + delta_z[i]

        self.channelUndulation = z_v
        self.thalweg = y - z_v + self.zd + zd


    def getThalweg(self):
        '''Return self.thalweg'''
        if self.thalweg is None:
            self.setThalweg()
        return self.thalweg


    def get_s_len(self):
        '''Return the length of centerline'''
        s_center = self.getCenterline_sn()
        return s_center[-1]


    def getDynamicPipeSlope(self):
        '''
        calculate the pipe slope of every point based on the sinuosity.
        '''
        x_v = self.getCenterline_x()
        s_v = self.getCenterline_sn()
        delta_x_v = np.concatenate((np.array([1]), x_v[1:] - x_v[0:-1]))
        delta_s_v = np.concatenate((np.array([1]), s_v[1:] - s_v[0:-1]))

        return self.x_slope * delta_x_v / delta_s_v

    # TODO: VALLEY CALC
    def getPipeSlope(self):
        '''Return pipe slope'''
        s_center = self.getCenterline_sn()
        if s_center[-1] == self.x_v[-1]:
            return self.x_slope
        return self.x_slope*self.x_len / self.get_s_len()


    def getSL(self):
        '''Return sinuosity of centerline'''
        return self.s_center[-1] / self.x_v[-1]


    def getHeight(self, direction):
        '''Return z value of most outer bank of river

        direction -- "left" or "right"
        '''
        if self.levels_z[direction] == []:
            return None
        return np.amax(self.levels_z[direction][-1])


    def setLevel(self, z_offset, z_start, y_offset, direction, yfun=None, innerPipe=False):
        '''Add one more level of level in a certain direction

        z_offset -- float; offset on z direction
        z_start -- array; z values of previous level
        y_offset -- float; offset on y direction
        direction -- "left" or "right"
        yfun -- function to calculate y values of level
        innerPipe -- boolean; if this is an inner pipes calculation
        '''
        z_offset = z_offset/self.dx
#        y_offset = y_offset/self.dx
        
        x_v = self.getCenterline_x()
        y_v = self.getCenterline_y()

        start = np.array([0]*len(x_v))
        prevLevel_x = x_v
        prevLevel_y = y_v

        level_x, level_y, level_n = self.addLevelOffset(x_v, y_v, start, y_offset, direction, self.getCenterline_sn(), yfun)
        innerPipePoints = {}
        level_x, level_y, ends_origin, innerPipePoints = self.levelCleanUp(level_x, level_y, x_v, y_v, direction, level_n, innerPipe)
        level_z = self.addTopOffset(level_x, level_y, prevLevel_x, prevLevel_y, z_start, z_offset)

        self.levels_n[direction].append(level_n)
        self.levels_x[direction].append(level_x)
        self.levels_y[direction].append(level_y)
        self.levels_z[direction].append(level_z)

        if innerPipe:
            if self.innerPipePoints == {}:
                self.innerPipePoints = innerPipePoints
            else:
                for x in innerPipePoints.keys():
                    if x in self.innerPipePoints:
                        self.innerPipePoints[x] = self.innerPipePoints[x].union(innerPipePoints[x])
                    else:
                        self.innerPipePoints[x] = innerPipePoints[x]


    def getLevel(self, dictkey, direction, ind=None):
        '''Return one level or list of levels in the given direction.
        
        dictkey -- 'x', 'y', 'z', or 'n'
        direction -- 'left' or 'right'
        ind -- index of the level wanted. If None, return all levels.
        '''
        levels = self.leveldict[dictkey][direction]
        if ind is None:
            return levels
        elif levels == []:
            return None
        else:
            if ind >= len(levels) or -1*ind > len(levels):
                print("Doesn't have", ind, "level in",direction,'in', dictkey)
                print('Return the last level instead.')
                return levels[-1]
            return levels[ind]

################################
    def createCenterline(self, fun):
        '''Calculate y values of centerline

        fun -- function used to calculate centerline
        '''
        x_center, y_center = fun(self.x_index)
        y_center = y_center/self.dx
        return x_center, y_center


    def calculate_s(self, x_v, y_v):
        '''Return s vector from xy to sn system.'''
        x1_v = x_v[:-1]
        x2_v = x_v[1:]
        y1_v = y_v[:-1]
        y2_v = y_v[1:]
        sqrt_v = np.vectorize(sqrt)
        ds_v = sqrt_v((y2_v-y1_v)**2 + (x2_v - x1_v)**2)
        s_v = [x_v[0]]
        si = 0

        for ds in np.nditer(ds_v, order="C"):
            si += ds
            s_v.append(si)
        return np.array(s_v)


    def addLevelOffset(self, x_v, y_v, start, minimum, direction, funInput, fun=None):
        '''Return x vector, y vector, and n vector for an additional level.

        x_v -- array; x values of the centerline
        y_v -- array; y values of the centerline
        start -- array; values of previous level
        minimum -- float; minimum width offset
        direction -- "left" or "right"
        funInput -- array; values used as input to later function
        fun -- function use to calculate the level

        Return:
        level_x -- array; x values of the new level.
        level_y -- array; y values of the new level.
        level_n -- array; n values of the new level.
        '''
        level_n = functions.offset_v(funInput, start, minimum, direction, fun)/self.dx
        level_x, level_y = functions.sn_to_xy_v(x_v, y_v, level_n)
        return level_x, level_y, level_n


    def levelCleanUp(self, level_x, level_y, previous_x, previous_y, direction, level_n, innerPipe):
        '''Clean up points of level cross over each other after it is transformed from sn to xy.

        level_x -- array; x values
        level_y -- array; y values
        direction -- "left" or "right"
        level_n -- array; n values

        Return:
        Clean up level_x and level_y.
        '''
        # Doesn't need to clean if no curvatue
        points_end = [(int(round(level_x[i])), int(round(level_y[i]))) for i in range(len(level_x))]

        x_max_current = np.amax(level_x)
        x_max_previous = np.amax(previous_x)
        x_max = ceil(max(x_max_current, x_max_previous))

        x_min_current = np.amin(level_x)
        x_min_previous = np.amin(previous_x)
        x_min = floor(min(x_min_current, x_min_previous))

        points_start = [(int(round(previous_x[i])), int(round(previous_y[i]))) for i in range(len(previous_x)) ]
        points_end_set = set(points_end)

        # compute covered area
        areaList_len = x_max - x_min + 3
        if x_min > 0:
            areaList_len = x_max + 3

        covered_area = self.calCoveredArea(areaList_len, points_start, points_end)

        innerPipePoints = {}
        if innerPipe:
            for x in range(x_min, x_max+1):
                innerPipePoints[x] = covered_area[x]

        out_x, out_y = self.retrieveOuterShape(covered_area, floor(x_min_current), ceil(x_max_current), direction, points_start, level_n, innerPipe, self.innerPipePoints)
        return out_x, out_y, points_end, innerPipePoints


    def buildVector(self, p1, p2):
        '''Return a list of (x, y) tuples that are crossed by line from point1 to point2

        Return:
        [(xi, yi)] -- xi, yi are rounded to 1
        '''
        vector = set()
        # if it is vertical line
        if p2[0] == p1[0]:
            ystart = min(p2[1], p1[1])
            yend = max(p2[1], p1[1])
            y = list(range(ystart, yend+1))
            vector = [(p1[0], yi) for yi in y]
        else:
        # if it is not a vertical line
            p1x, p1y = p1[0], p1[1]
            p2x, p2y = p2[0], p2[1]

            s = functions.slope(p1x, p1y, p2x, p2y)
            step = ceil(functions.pointDist(p1, p2))
            dx = (p2x-p1x)/step

            for t in range(step+1):
                dot_x = p2x - t*dx
                dot_y = p2y - s*t*dx

                vector.add((round(dot_x), round(dot_y)))
        vector = list(vector)
        vector.sort()

        return vector


    def addTopOffset(self, newLevel_x, newLevel_y, prevLevel_x, prevLevel_y, prevLevel_z, offset):
        ''' Calculate the height of a new added level.

            The height is calculated from its previous level and
                the slope of the pipe.
        '''
        newLevel_z = np.empty(len(newLevel_x))

        for i in range(len(newLevel_x)):
            x = newLevel_x[i]
            y = newLevel_y[i]

            dist = np.square(prevLevel_x - x) + np.square(prevLevel_y - y)
            minIndex = np.argmin(dist)
            newLevel_z[i] = prevLevel_z[minIndex] + offset

        return newLevel_z


    def commonX(self, v0, v1, covered):
        xi0 = set([p[0] for p in v0])
        xi1 = set([p[0] for p in v1])
        x_inter = xi0.intersection(xi1)

        for x in x_inter:
            y0 = [p[1] for p in v0 if p[0] == x]
            y1 = [p[1] for p in v1 if p[0] == x]
            y = y0 + y1
            ymin = min(y)
            ymax = max(y)

            for y in range(ymin, ymax+1):
                covered[x].add(y)
        return covered


    def calCoveredArea(self, areaListLength, points_start, points_end):
        '''
        Calculate area covered by a list of vectors.
        '''
        covered_area = [set() for i in range(areaListLength)]
        p1_pre = points_start[0]
        p2_pre = points_end[0]
        v0 = self.buildVector(p1_pre, p2_pre)
        for (x, y) in v0:
            covered_area[x].add(y)

        for i in range(1, len(points_start)):
            p1 = points_start[i]
            p2 = points_end[i]

            v0 = self.buildVector(p1_pre, p2_pre)
            v1 = self.buildVector(p1, p2)
            for (x, y) in v1:
                covered_area[x] = covered_area[x]
                covered_area[x].add(y)

            v2 = self.buildVector(p2_pre, p2)
            for (x, y) in v2:
                covered_area[x] = covered_area[x]
                covered_area[x].add(y)

            covered_area = self.commonX(v0, v1, covered_area)
            covered_area = self.commonX(v1, v2, covered_area)
            covered_area = self.commonX(v0, v2, covered_area)

            p1_pre = p1
            p2_pre = p2

        for i in range(1, len(covered_area)-1):
            ys = covered_area[i]
            pre_ys = covered_area[i-1]
            post_ys = covered_area[i+1]

            inter_ys = pre_ys.intersection(post_ys)
            for t in inter_ys:
                covered_area[i].add(t)

            ys = list(ys)
            ys.sort()

            if len(ys) > 1:
                for t in range(1, len(ys)):
                    y_pre = ys[t-1]
                    y = ys[t]
                    if y - y_pre == 2:
                        covered_area[i].add(y-1)

        return covered_area

    
    def extractInnerPipePoints(self, innerPipePoints):
        '''Convert points in a dictionary to a set with points in (x, y) format.

        innerPipePoints - dictionary;
                            key: int; x value
                            value: list; y values

        Return:
        innerPoints - set; element: (x, y)
        '''
        innerPoints = set()
        for x in innerPipePoints.keys():
            ys = innerPipePoints[x]
            for y in ys:
                innerPoints.add((x, y))
        return innerPoints


    def checkHighCurvArea(self, x_min, x_max, points_start, level_n):
        '''Check if there are highcurve area and return a boolean list to indicate that.

        x_min - int; x value we start from to check
        x_max - int; x value we end on to check
        points_start - set; points in (x, y) format we want to check
        level_n - array; values of n values of level

        return:
        highcurvs - list of boolean; True means at this point there is a highcurv.
        '''
        x_min_0 = min(x_min, 0)
        highcurvs = [set() for i in range(x_max - x_min_0 + 1)]
        
        for (x, y) in points_start:
            if x < x_min or x > x_max:
                continue
            highcurvs[x].add(y)

        highcurvs = [list(ys) for ys in highcurvs]
        for i in range(len(highcurvs)):
            yslist = highcurvs[i]
            highcurvs[i] = False
            if len(yslist) <= 1:
                continue

            yslist.sort()
            for t in range(1, len(yslist)):
                if yslist[t] - yslist[t-1] > 4:
                    highcurvs[i] = True
                    break
        
        average_n = int(np.average(level_n))

        if True in highcurvs:
            for i in range(len(highcurvs)):
                checkTrue = highcurvs[i]
                if not checkTrue:
                    for t in range(max(x_min, i-average_n), min(x_max, i+average_n)+1):
                        if highcurvs[t]:
                            highcurvs[i] = True
                            break
        return highcurvs


    def retrieveOuterShape(self, covered_area, x_min, x_max, direction, points_start, level_n, innerPipe, innerPipePoints):
        '''
        Retrieve x values and y values for the outer shape of the covered_area

        covered_area - list of lists; index is x value, lists are y values
        x_min - int; x value we start from to draw the outer shape
        x_max - int; x value we end on to draw the outer shape
        direction - 'left' or 'right'
        points_start - set of tuples; contains starting points as in (x, y) format
        level_n - array; values of n values of level
        innerPipePoints - dictionary;
                            key: int; x value
                            value: list; y values

        Return:
        out_x - array; x values of level
        out_y - array; y values of level
        '''
        # extract points from innerPoints dictionary
        innerPoints = self.extractInnerPipePoints(innerPipePoints)

        # check if there are highcurv area
        # highcurv area means one x value with multiple y values
        highcurvs = self.checkHighCurvArea(x_min, x_max, points_start, level_n)

        # start retrieving outer shape.
        out_x = []
        out_y = []
        startSet = set(points_start)

        for x in range(x_min, x_max+1):
            # if it is empty set
            if not covered_area[x]:
                continue
            # if it is not empty set
            yli = list(covered_area[x])
            if direction == 'left':
                yli.sort(reverse=True)
            else:
                yli.sort()
            if (x, yli[0]) not in innerPoints and ((x, yli[0]) not in startSet):
                out_x.append(x)
                out_y.append(yli[0])

            if len(yli) < 3:
                continue
            y_pre = yli[0]
            for i in range(1, len(yli)-1):

                y = yli[i]
                if (x, y) in innerPoints or (x, y) in startSet:
                    y_pre = y
                    continue

                if innerPipe:
                    if ((x, y+1) in startSet or (x, y+2) in startSet
                            or (x, y-1) in startSet or (x, y-2) in startSet):
                        y_pre = y
                        continue

                    if y in covered_area[x-1] and y in covered_area[x-2]:
                        if y not in covered_area[x+1] and y not in covered_area[x+2]:
                            if (x-1, y) not in startSet and (x-2, y) not in startSet:
                                out_x.append(x)
                                out_y.append(y)
                                y_pre = y
                                continue
                    elif y not in covered_area[x-1] and y not in covered_area[x-2]:
                        if y in covered_area[x+1] and y in covered_area[x+2]:
                            if (x+1, y) not in startSet and (x+2, y) not in startSet:
                                out_x.append(x)
                                out_y.append(y)
                                y_pre = y
                                continue
                    else:
                        y_pre = y
                        continue

                if abs(y_pre-y) > 2:
                    out_x.append(x)
                    out_y.append(y)
                    if i != 1 and (x, y_pre) not in startSet and ((x, y_pre) not in innerPoints):
                            out_x.append(x)
                            out_y.append(y_pre)
                y_pre = y
            
            # only high curvature lines need to check end point
            if highcurvs[x]:
                if ((x, yli[-1]) not in startSet
                        and (x-1, yli[-1]) not in startSet
                        and (x+1, yli[-1]) not in startSet):
                    out_x.append(x)
                    out_y.append(yli[-1])
                

        out_x = np.array(out_x)
        out_y = np.array(out_y)

        return out_x, out_y
