'''This module simulate a river-valley system based on user inputs.

Usage:
python3 riverbuilder <path.to.input.txt> <outputFolderName>

path.to.input.txt -- Absolute or relative path to an input file that contains
                    all parameters needed to build a river.
outputFolderName -- Name of the folder that output files will be stored in. If
                    the folder doesn't exist, then it will be created.

Its overall mechanism is:
    1. Parse inputs from input file.
    2. Check and convert inputs.
    3. Build a corresponding river.
    4. Build a corresponding valley.
'''

from .functions import *
from .cChannel import Channel
from .cValley import Valley
from math import pi
from decimal import Decimal
import sys
import csv
import os
import re
import copy
import traceback
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime


ALLOWFUN = set(["MASK", "SIN", "COS", "LINE", "PERLIN", "SINSQ", "COSSQ", "CNOIDAL", "STEP", "HIGHCURV"])
FUNPARANUM = {"SIN":4, "COS":4, 'SINSQ':4, 'COSSQ':4, "LINE":3, "PERLIN":4, "CNOIDAL":5, "STEP":5,
        "HIGHCURV": 5}
XSHAPES = set(['AU', 'SU', 'EN', 'CF', 'PY', 'AF', 'DT', 'TU'])
ADDON = {'BEG':[], 'CD':[]}


def defaultFunction(x):
    '''The default function for any function calculations, which returns a line without fluctuation.'''
    return x, np.array([0]*len(x))


def readTXTFunctionFile(nameId, paraDict, paraName, textName):
    try:
        f = open(textName, "r")
    except IOError:
        log = "Error! "+textName+" is not found in current folder:\n"+os.getcwd()+"\nProgram exits.\n"
        print(log)
        sys.exit()

    lines = f.readlines()
    funID = 0
    funList = []
    for line in lines:
        if line.startswith('#') or '=' not in line:
            continue
        [name, val] = line.strip().split('=')

        if name in paraDict:
            name = re.split('\d+', name)[0]
            name = name + '0' + str(nameId) + str(funID)
            funID += 1

        paraDict[name] = val
        if not name.startswith('MASK'):
            funList.append(name)

    funList = '+'.join(funList)

    if paraName in paraDict:
        paraDict[paraName] = paraDict[paraName] + '+' + funList
    else:
        paraDict[paraName] = funList

    return paraDict


def fileParser(fname):
    '''Parse a file, return a dictionary

    fname -- str: path to file

    Return:
    outdict -- dictionary, both key and value are string

    Exception:
    IOError -- program will exit.
    '''
    outdict = {}
    try:
        f = open(fname, "r")
    except IOError:
        log = "Error! "+fname+" is not found in current folder:\n"+os.getcwd()+"\nProgram exits.\n"
        print(log)
        sys.exit()

    lines = f.readlines()
    addonSet = ADDON.keys()
    p = re.compile('[A-Z]+')
    nameID = 0
    for line in lines:
        if line.startswith('#') or '=' not in line:
            continue
        [name, val] = line.strip().split('=')
        # check if it is an add-on:
        m = p.match(name)
        if m.group() is not None and m.group() in addonSet:
            ADDON[m.group()].append(val)
            continue

        if val.endswith('txt'):
            outdict = readTXTFunctionFile(nameID, outdict, name, val)
            nameID += 1
            continue

        if name not in outdict:
            outdict[name] = val
        else:
            outdict[name] = outdict[name]+"+"+val

    return outdict


#def fileParser(fname):
#    '''Parse a file, return a dictionary
#
#    fname -- str: path to file
#
#    Return:
#    outdict -- dictionary, both key and value are string
#
#    Exception:
#    IOError -- program will exit.
#    '''
#    outdict = {}
#    try:
#        f = open(fname, "r")
#    except IOError:
#        log = "Error! "+fname+" is not found in current folder:\n"+os.getcwd()+"\nProgram exits.\n"
#        print(log)
#        sys.exit()
#
#    lines = f.readlines()
#    addonSet = ADDON.keys()
#    p = re.compile('[A-Z]+')
#    for line in lines:
#        if line.startswith('#') or '=' not in line:
#            continue
#        [name, val] = line.strip().split('=')
#        # check if it is an add-on:
#        m = p.match(name)
#        if m.group() is not None and m.group() in addonSet:
#            ADDON[m.group()].append(val)
#            continue
#
#        if name not in outdict:
#            outdict[name] = val
#        else:
#            outdict[name] = outdict[name]+"+"+val
#
#    return outdict
#
#
def paraCheck(fdict, name, defaultVal, valType, sign=0):
    '''Check if the value for a key in a dictionary is valid, return updated dictionary

    fdict -- dictionary needs to be checked and updated
    name -- the key in fdict need to be checked
    defaultVal -- default value for name in fdict when checked failed
    valType -- str that indicate appropriate type of value for name in fdict
    sign -- -1, 0, or 1 represents whether value should be positive, real, or negative

    Return:
    fdict -- updated dictionary with correct value for name
    log -- str that is empty or alert information
    '''
    log = ""
    if name not in fdict:
        fdict[name] = defaultVal
        log += 'Alert! '+name+' is not defined by user, use default value '+str(defaultVal)+' instead.\n'
        return fdict, log

    if valType == 'str':
        return fdict, log

    try:
        val = float(fdict[name])
    except:
        log += "Alert! Can't not convert value for "+name+" to a number.\n Use default value "+defaultVal+" instead.\n"
        fdict[name] = defaultVal
    if valType == "int":
        fdict[name] = int(fdict[name])
    elif valType == "float":
        fdict[name] = float(fdict[name])

    if sign != 0 and fdict[name]*sign < 0:
        log += "Alert! The sign of value of "+name+" is incorrect. Change to the opposite.\n"
        fdict[name] = -1* fdict[name]

    return fdict, log


def getfunParas(val, funDict):
    '''Parse parameter set string for a function, return a correct callable function.

    val -- parameter set string

    Constrains:
    This function should only be called when both name and val are in correct form.
    '''
    val = val[1:-1].split(",")
    mask = val.pop().strip()
    if mask.startswith('MASK'):
        if mask in funDict:
            mask = funDict[mask]
        else:
            mask = ['ALL']

    for i in range(len(val)):
        num = val[i]
        num = num.strip()
        exp = re.split(r'[/*+-]', num)
        if exp[0].strip() == '':
            exp[0] = 0
        exp = [pi if x == "pi" else float(x) for x in exp]
        if len(exp) == 1:
            val[i] = float(exp[0])
            continue
        if '+' in num:
            val[i] = exp[0] + exp[1]
        elif '-' in num:
            val[i] = exp[0] - exp[1]
        elif '*' in num:
            val[i] = exp[0] * exp[1]
        elif '/' in num:
            val[i] = exp[0] / exp[1]

    val.append(mask)
    return val


def forgeFun(name, val):
    '''
    name -- name of function
    '''
    if name == "SIN":
        return lambda x: sin_v(x, *val)
    elif name == "COS":
        return lambda x: cos_v(x, *val)
    elif name == "SINSQ":
        return lambda x: sin_sq(x, *val)
    elif name == "COSSQ":
        return lambda x: cos_sq(x, *val)
    elif name == "LINE":
        return lambda x: line_v(x, *val)
    elif name == "PERLIN":
        return lambda x: perlin(x, *val)
    elif name == "CNOIDAL":
        return lambda x: cnoidal(x, *val)
    elif name == "STEP":
        return lambda x: step(x, *val)
    elif name == "HIGHCURV":
        return lambda x: highCurv(x, *val)


def buildFunDict(fdict):
    '''Given a dictionary, extract and parse all allowed functions in it and return a new dictionary

    fdict -- dictionary that may contains allowed functions which are in str type.

    Return:
    funDict -- dictionary that contains allowed functions which values are callable functions.
    log -- str that record all successfully parsed functions and alert msg.

    Constrains:
    Need to defined allowFun and modify funValCheck and forgeFun when the code is changed.
    '''
    log = ""
    p = re.compile('[A-Z]+')
    funDict = {}

    notMaskFuns = []
    for name, val in fdict.items():
        if not isinstance(val, str):
            continue
        name = name.strip()
        val = val.strip()
        m = p.match(name)
        if (val.startswith("(") and
                  val.endswith(")") and
                  m is not None and
                  m.group() in ALLOWFUN):

            if m.group() == 'MASK':
                funDict[name] = val[1:-1].split(',')
                continue

            if not funValCheck(m.group(), val):
                log += "Can't parse function "+ name +'\n'
                continue

            notMaskFuns.append([name, m.group(), val])

    for [name, funType, val] in notMaskFuns:
        funParas = getfunParas(val, funDict)
        funDict[name] = forgeFun(funType, funParas)

    log += "User defined funtions are:\n"
    for fun in funDict:
        log += fun+' '+fdict[fun]+'\n'

    return funDict, log


def funValCheck(name, val):
    '''Check if parameter set string is correct for a function, return boolean

    name -- name of function
    val -- parameter set string in form (n1, n2, ...)
    '''
    val = val[1:-1].split(",")

    if len(val) != FUNPARANUM[name]:
        return False

    val[-1] = val[-1].strip()
    if val[-1].startswith("MASK"):
        val.pop()

    for num in val:
        num = num.strip()
        p = re.compile("(\d+.?\d*|pi)")
        splitlist = re.split(r'[/*+-]', num)
        for i in range(len(splitlist)):
            part = splitlist[i]
            part = part.strip()
            if i == 0 and part == '':
                continue
            m = p.match(part)
            if not m:
                return False
            if m.span() != (0, len(part)):
                return False
    return True


def printPara(name, fdict):
    '''Return a name value pair in dictionary in string form.'''
    return name+' is set to: '+str(fdict[name])+'\n'


def inputCheck(fdict):
    '''Check keys and values in dictionary, return updated dictionaries that
    contains all important information to build river and valley.

    fdict -- dictionary need to be checked

    Return:
    fdict -- updated dictionary with correct key and value pairs
    funDict -- dictionary with values as callable functions
    log -- string format log
    '''
    log = ''
    fdict, info = paraCheck(fdict, "Datum", 10, "float")
    log += info
    log += printPara("Datum", fdict)
    fdict, info = paraCheck(fdict, "Length", 1000, "float", 1)
    log += info
    log += printPara("Length", fdict)
    fdict, info = paraCheck(fdict, "X Resolution", 1, "float", 1)
    log += info
    log += printPara("X Resolution", fdict)
    fdict, info = paraCheck(fdict, "Channel XS Points", 21, "int", 1)
    log += info
    log += printPara("Channel XS Points", fdict)
    fdict, info = paraCheck(fdict, "Valley Slope (Sv)", 0.001, "float")
    log += info
    log += printPara("Valley Slope (Sv)", fdict)
    fdict, info = paraCheck(fdict, "Critical Shields Stress (t*50)", 0.06, "float")
    log += info
    log += printPara("Critical Shields Stress (t*50)", fdict)
    fdict, info = paraCheck(fdict, "Inner Channel Lateral Offset Minimum", 10, "float", 1)
    log += info
    log += printPara("Inner Channel Lateral Offset Minimum", fdict)
    fdict, info = paraCheck(fdict, "Inner Channel Depth Minimum", 0, "float", 0)
    log += info
    log += printPara("Inner Channel Depth Minimum", fdict)
    fdict, info = paraCheck(fdict, "Median Sediment Size (D50)", 0.01, "float", 1)
    log += info
    log += printPara("Median Sediment Size (D50)", fdict)
    fdict, info= paraCheck(fdict, "Cross-Sectional Shape", "SU", "str")
    log+= info
    log+=printPara("Cross-Sectional Shape", fdict)
    fdict, info = paraCheck(fdict, "Left Valley Boundary Lateral Offset Minimum", 10, "float", 1)
    log += info
    log += printPara("Left Valley Boundary Lateral Offset Minimum", fdict)
    fdict, info = paraCheck(fdict, "Left Valley Boundary Height Offset", 20, "float", 1)
    log += info
    log += printPara("Left Valley Boundary Height Offset", fdict)
    fdict, info = paraCheck(fdict, "Right Valley Boundary Lateral Offset Minimum", 10, "float", 1)
    log += info
    log += printPara("Right Valley Boundary Lateral Offset Minimum", fdict)
    fdict, info = paraCheck(fdict, "Right Valley Boundary Height Offset", 20, "float", 1)
    log += info
    log += printPara("Right Valley Boundary Height Offset", fdict)
    fdict, info = paraCheck(fdict, "Inner Channel Average Bankfull Width", None, "float", 1)
    log += info
    log += printPara("Inner Channel Average Bankfull Width", fdict)
    fdict, info = paraCheck(fdict, "Inner Channel Average Bankfull Depth", None, "float", 1)
    log += info
    log += printPara("Inner Channel Average Bankfull Depth", fdict)
    fdict, info = paraCheck(fdict, "River Slope", None, "float", 1)
    log += info
    log += printPara("River Slope", fdict)
    fdict, info = paraCheck(fdict, "Smooth", 0, "int", 1)
    log += info
    log += printPara("River Slope", fdict)
    log += ''

    funDict, info = buildFunDict(fdict)
    log += info
    return fdict, funDict, log


def funParser(funString, funDict, defaultFun):
    '''Parse a string of different functions, return a callable aggregate function.

    funString -- string in format "SIN1+COS2+LINE3"
    funDict -- dictionary with key as individual function name, 
                and value as corresponding callable function

    Return:
    outFun -- callable aggregated function
    log -- string log
    '''
    log = ''
    funs = funString.split("+")

    notIn = True
    highcurvFlag = False
    for i in range(len(funs)):
        fun = funs[i]
        if fun in funDict:
            notIn = False
        else:
            log += "Alert! Can't find function "+fun+" in user-defined functions. Ignore function "+fun+'.\n'
        #switch the high curv functions to the top
        if fun.startswith('HIGHCURV'):
            funs.pop(i)
            if highcurvFlag:
                log += "Error! Can only have one HIGHCURV function! Extra ones will be ignored!"
                continue
            funs.insert(0, fun)
            highcurvFlag = True

    if notIn:
        return defaultFun, log

    def outFun(x):
        outSum = 0
        for fun in funs:
            if fun in funDict:
                x, out = funDict[fun](x)
                outSum += out

        return x, outSum

    return outFun, log


def buildFun(name, fdict, funDict, fun):
    '''Return an appropriate callable function for a given name.

    name -- the thing that we want to assign a function
    fdict -- dictionary to look up what kind of functions are needed for name
    funDict -- dictionary with key as individual function name,
                and value as corresponding callable function
    fun -- default function if the check is failed

    Return:
    outFun -- callable function for name
    log -- string log
    '''
    log = ''
    outFun = fun
    if name in fdict:
        outFun, info = funParser(fdict[name], funDict, outFun)
        log += info
    else:
        log += "Alert! Can't find function "+name+" in user-defined functions. Ignore function "+name+'.\n'
    return outFun, log


def addLevels(pattern, fdict, funDict, default_fun,direction, obj):
    '''Add levels to river or valley in a given direction.

    pattern -- regex pattern that match the name of levels
    fdict -- dictionary that gives information of levels wanted to add
    funDict -- dictionary contain individaul function information
    default_fun -- default function if check or parse fails
    direction -- string "left" or "right"
    obj -- a Channel object or Valley object

    Return:
    obj -- an updated obj with levels added
    log -- string log
    '''
    log = ''
    levels = []

    for name, val in fdict.items():
        m = pattern.match(name)
        if m:
            levels.append(name)

    levels.sort()

    for name in levels:
        funName = name[:-22] + "Function"
        fun, info = buildFun(funName, fdict, funDict, default_fun)
        log += info

        fdict, info = paraCheck(fdict, name, 10, "float", 1)
        log += info
        y_offset = fdict[name]

        hightName = name[:-22]+"Height Offset"
        fdict, info = paraCheck(fdict, hightName, 10, "float")
        log += info
        z_offset = fdict[hightName]

        thalweg = obj.getThalweg()
        if hasattr(obj, 'channel'):
            z_pre = thalweg
        else:
            orig_thalweg = thalweg+obj.channelUndulation
            thalweg_max = np.amax(orig_thalweg)
            z_pre = thalweg_max - obj.channelUndulation

        obj.setLevel(z_offset, z_pre, y_offset, direction, fun)
        if funName in fdict:
            log += "Creating "+name[:-22]+"with function: "+str(fdict[funName])+'\n'
        else:
            log += "Creating "+name[:-22]+"with constant width: "+str(fdict[name])+'\n'

    return obj, log


def loopParameter(para, target, obj, buildfun, calfun):
    '''Rebuild the channel again and again until a parameter reach target value.
    
    para - parameter that will be modified every loop
    target - target value of the parameter
    minimum - minimum offset set for the parameter
    obj - channel or river
    buildfun - function that used to build the parameter, ex: channel.createInnerChannel
    calfun - function that used to calculate the parameter

    Return:
    obj - object with a target paraName
    log - string log
    '''
    log = ''
    count = 0
    decNum = abs(Decimal(str(target)).as_tuple().exponent)

    objTemp = copy.deepcopy(obj)
    objTemp = buildfun(objTemp, para)
    out = round(calfun(objTemp), decNum)
    sign = target - out
    increment = target - out
    flag = False

    while out != target and count < 100:
        para += increment

        if para <= 0 and flag:
            return -1, log
        elif para <= 0 and not flag:
            para = 0
            flag = True
        else:
            flag = False

        objTemp = copy.deepcopy(obj)
        objTemp = buildfun(objTemp, para)
        out = round(calfun(objTemp), decNum)

        if sign * (target - out) < 0:
            para -= increment
            increment = increment/2
        else:
            increment = target - out

        count += 1

    if out != target:
        log = 'Looping reaches a limit of 100 times. Stop.\n'
    else:
        log = 'Loop '+str(count+1)+' times to get the ideal value.\n'
    print(log)

    return para, ''


def buildChannel(fdict, funDict):
    '''Create a channel basing on information given in fdict and funDict.

    fdict -- dictionary contains parameter information for a channel
    funDict -- dictionary contains callable user defined functions

    Return:
    c -- channel that been built
    log -- string log
    '''
    log = ''
    c = Channel(fdict["Length"], fdict["Inner Channel Lateral Offset Minimum"],\
            fdict["Valley Slope (Sv)"], fdict["X Resolution"], fdict["Datum"])
    nPoints = fdict['Channel XS Points']
    c.setXShapePoints(nPoints)


    valleyfun, info = buildFun("Valley Centerline Function", fdict, funDict, defaultFunction)
    reshape=True

    fun, info = buildFun("Meandering Centerline Function", fdict, funDict, defaultFunction)
    log += info

    if valleyfun == defaultFunction:
        reshape = False
        log += 'Reshape not needed for river centerline.\n'

    if fdict['River Slope'] is not None:
        if fdict['River Slope'] >= fdict['Valley Slope (Sv)']:
            c.setCenterline(fun)
            log += 'Error! River Slope can not be bigger than Valley Slope!\n'
            print('Error! River Slope can not be bigger than Valley Slope!')
        else:
            print('')
            print('Start looping to get target river slope...')

            para = 1

            rslope = fdict['River Slope']
            count = 0
            decNum = abs(Decimal(str(rslope)).as_tuple().exponent)

            channelTemp = copy.deepcopy(c)
            
            def centerlinefun(x):
                x, y = fun(x)
                return x, para*y

            channelTemp.setCenterline(centerlinefun)

            if reshape:
                print('reshape: Yes')
                channelTemp.shapeCenterline(valleyfun)

            out = channelTemp.getPipeSlope()
            increment = rslope/out

            while out != rslope and count < 100:
                para /= increment

                channelTemp = copy.deepcopy(c)
                channelTemp.setCenterline(centerlinefun)

                if reshape:
                    channelTemp.shapeCenterline(valleyfun)

                out = round(channelTemp.getPipeSlope(), decNum)

                increment = rslope/out
                count += 1

            if out != rslope:
                print('Looping reaches a limit of 100 times. Stop.')
            else:
                print('Loop '+str(count+1)+' times to get the ideal value.')

            c.setCenterline(centerlinefun)
    else:
        c.setCenterline(fun)

    if reshape:
        c.shapeCenterline(valleyfun)

    c.smoothCenterline(fdict['Smooth'])
    log += "Creating Meandering Center line with Function:"+fdict.get("Meandering Centerline Function", "None")+'\n'

    if fdict["Inner Channel Depth Minimum"] != 0:
        c.setHbfManual(fdict["Inner Channel Depth Minimum"])
        log += "Use user defined Inner Channel Depth Minimum.\n"

    fun, info = buildFun("Centerline Curvature Function", fdict, funDict, None)
    log += info
    if fun:
        c.setCurvature(fun)
        log += "Use user defined Centerline Curvature Function:"+fdict["Centerline Curvature Function"]+'\n'

    leftfun, info = buildFun("Left Inner Bank Function", fdict, funDict, None)
    log += info
    rightfun, info = buildFun("Right Inner Bank Function", fdict, funDict, None)
    log += info
    thalfun, info = buildFun("Thalweg Elevation Function", fdict, funDict, defaultFunction)
    log += info
    print('')

    if fdict['Inner Channel Average Bankfull Width'] is not None:
        print('Start looping to get target Inner Channel Average Bankfull Width...')

        def loopFun(channel, para):
            channel.wbf_min = para
            channel.createInnerChannel(leftfun, rightfun, thalfun)
            return channel

        def calFun(channel):
            return channel.getAveWbf()

        para = fdict["Inner Channel Lateral Offset Minimum"]
        para, info = loopParameter(para, fdict['Inner Channel Average Bankfull Width'], c, loopFun, calFun)
        if para == -1:
            log += 'Cannot reach target Inner Channel Average Bankfull Width with current function settings. Please modify the functions.'
        else:
            log += info
            c.wbf_min = para

    if fdict['Inner Channel Average Bankfull Depth'] is not None:
        print('Start looping to get target Inner Channel Average Bankfull Depth...')

        def loopFun(channel, para):
            channel.hbf = para
            channel.createInnerChannel(leftfun, rightfun, thalfun)
            return channel

        def calFun(channel):
            return channel.getAveHbf()

        para = fdict["Inner Channel Depth Minimum"]
        para, info = loopParameter(para, fdict['Inner Channel Average Bankfull Depth'], c, loopFun, calFun)
        if para == -1:
            log += 'Cannot reach target Inner Channel Average Bankfull Depth with current function settings. Please modify the functions.'
        else:
            log += info
            c.hbf = para

    c.createInnerChannel(leftfun, rightfun,  thalfun)
    log += "Creating Inner Channel Banks with left bank function: "+fdict.get("Left Inner Bank Function", "None")+'\n'
    log += "                             with right bank function: "+fdict.get("Right Inner Bank Function", "None")+'\n'
    log += "                             with thalweg elevation function: "+fdict.get("Thalweg Elevation Function", "None")+'\n'

    pattern = re.compile('L[\d]+ Outer Bank Lateral Offset Minimum')
    c, info = addLevels(pattern, fdict, funDict, None, 'left', c)
    log += info

    pattern = re.compile('R[\d]+ Outer Bank Lateral Offset Minimum')
    c, info = addLevels(pattern, fdict, funDict, None, 'right', c)
    log += info

    #################### Cross-Sectional Shape Calculation ################################
    ckey = 'Cross-Sectional Shape'
    if ckey not in fdict:
        log += 'Alert! Cross-Sectional Shape not specified! Use asymmetric shape as default.\n'
        fdict[ckey] = 'AU'
    
    if fdict[ckey] not in XSHAPES:
        log += 'Alert! Cross-Sectional Shape value not recognizable! User input: '+ str(fdict[ckey])+'\nUse asymmetric shape as default.\n'
        fdict[ckey] = 'AU'

    if fdict[ckey] == 'AU':
        c.setXShape()
    elif fdict[ckey] == 'SU':
        copyCurvature = copy.copy(c.getDynamicCurv())
        c.dynamicCurv = c.dynamicCurv*0
        c.setXShape()
    elif fdict[ckey] == 'CF':
        CF_a = fdict.get("CF_a")
        CF_b = fdict.get("CF_b")
        CF_c = fdict.get("CF_c")
        c.setXShapeCF(CF_a,CF_b,CF_c)
        c.cross_section_type = 'CF'
    elif fdict[ckey] == 'PY':
        c.setXShapePY()
        c.cross_section_type = 'PY'
    elif fdict[ckey] == 'AF':
        d1 = fdict.get("Depth1")
        d2 = fdict.get("Depth2")
        ang1 = fdict.get("Angle1")
        ang2 = fdict.get("Angle2")
        c.setXShapeAF(d1, d2, ang1, ang2)
        c.cross_section_type = 'AF'
    elif fdict[ckey] == 'DT':
        c.setXShapeDT()
        c.cross_section_type = 'DT'
    elif fdict[ckey] == 'TU':
        c.setXShapeTU()
        c.cross_section_type = 'TU'
    else:
        fdict, info = paraCheck(fdict, "TZ(n)", 1, "int", 1)
        log += info
        log += printPara("TZ(n)", fdict)
        if fdict['TZ(n)'] > nPoints:
            log += 'Alert! TZ(n) value is not valid, set to Channel XS Points.\n'
            fdict['TZ(n)'] = nPoints
        c.setXShape(fdict['TZ(n)'])
        c.setTZ(fdict['TZ(n)'])

    #################### Bed Roughness ################################
    if 'PBR' in fdict:
        log += 'Adding perlin roughness to bed.\n'
        fdict, info = paraCheck(fdict, "PBR", 5, "float")
        c.perlinThalweg(fdict['PBR'])

    return c, log


def addBedElement(c, paras):
    log = ''
    paras = paras[1:-1]
    paras = paras.split(',')
    paras = [i.strip() for i in paras]

    if len(paras) != 4:
        log += "Number of parameters given to BEG doesn't equal to 4.\n"
        log += 'Given parameters should be in form: (val1, val2, val3, val4)\n'
        return log

    [num, size_mean, size_std, height] = paras
    try:
        num = int(num)
        size_mean = int(size_mean)
        size_std = int(size_std)
        height = float(height)
    except ValueError:
        log += "Cannot parsed given parameters in BEG.\n"
        return log

    c.addBoulders(num, size_mean, size_std, height)
    return log


def addCheckDam(c, paras):
    log = ''
    paras = paras[1:-1]
    paras = paras.split(',')
    paras = [i.strip() for i in paras]

    if len(paras) != 3:
        log += "Number of parameters given to CD doesn't equal to 3.\n"
        log += 'Given parameters should be in form: (val1, val2, val3)\n'
        return log

    [loc, height, thickness] = paras
    try:
        loc = int(loc)
        height = int(height)
        thickness = int(thickness)
    except ValueError:
        log += "Cannot parsed given parameters in BEG.\n"
        return log

    c.addCheckDam(loc, height, thickness)
    return log

def addChannelElements(c, addon):
    ADDONMETHOD = {'BEG': addBedElement,
                    'CD': addCheckDam}

    log = ''
    for key, values in ADDON.items():
        if values != []:
            for value in values:
                log += ADDONMETHOD[key](c, value)
    return c, log


def buildValley(fdict, funDict, channel):
    '''Create a valley basing on information given in fdict and funDict.

    fdict -- dictionary contains parameter information for a valley
    funDict -- dictionary contains callable user defined functions
    channel -- channel that will be passed to valley

    Return:
    valley -- valley that been built
    log -- string log
    '''
    log = ''
    valley = Valley(channel)

    fun, info = buildFun("Valley Centerline Function", fdict, funDict, defaultFunction)
    log += info
    valley.setCenterline(fun)
    pattern = re.compile('(R[\d]+ Valley Breakline Lateral Offset Minimum)')
    valley, info = addLevels(pattern, fdict, funDict, None, 'right', valley)
    log += info

    pattern = re.compile('(L[\d]+ Valley Breakline Lateral Offset Minimum)')
    valley, info = addLevels(pattern, fdict, funDict, None, 'left', valley)
    log += info

    lboffset = fdict['Left Valley Boundary Lateral Offset Minimum']
    lbheight = fdict['Left Valley Boundary Height Offset']
    valley.setValleyBoundary(lbheight, lboffset, 'left', None)

    rboffset = fdict['Right Valley Boundary Lateral Offset Minimum']
    rbheight = fdict['Right Valley Boundary Height Offset']
    valley.setValleyBoundary(rbheight, rboffset, 'right', None)
    return valley, log


def plotLevels(ax, xdict, ydict, dx, labelend, col):
    '''Plot levels to a figure.

    ax - the ax to draw to 
    xdict - dictionary of values of x values of levels
    ydict - dictionary of values of y values of levels
    labelend - label added to each level, a string
    col - col of dots
    '''
    for i in range(len(xdict['left'])):
        x = xdict['left'][i]*dx
        y = ydict['left'][i]*dx
        if labelend == 'V':
            ax.scatter(x, y, c='C'+str(col+i), marker='.', label='L'+labelend+str(i+1))
        else:
            ax.scatter(x, y, c='C'+str(col+i), marker='.', label='L'+labelend+str(i))

    for i in range(len(xdict['right'])):
        x = xdict['right'][i]*dx
        y = ydict['right'][i]*dx
        if labelend == 'V':
            ax.scatter(x, y, c='C'+str(col+i), marker='_', label='R'+labelend+str(i+1))
        else:
            ax.scatter(x, y, c='C'+str(col+i), marker='_', label='R'+labelend+str(i))


###########################################################
def buildRiver(fname, outfolder, log):
    '''
    It parse parameters in inputfile, then output to outputfolder.

    fname - inputfile name, end with '.txt'.
    outfolder - output folder name.
    log - additional information want to added to log file.
    '''
    try:
        print("Start Parsing Inputs...")
        log += '\n'
        paraDict = fileParser(fname)

        paraDict, funDict, info = inputCheck(paraDict)
        log += info
        log += '\n'

        print("Start Creating River Channel...")
        t1 = datetime.now()
        channel, info = buildChannel(paraDict, funDict)
        log += info
        log += '\n'

        print("Start adding add-ons to Channel...")
        channel, info = addChannelElements(channel, ADDON)
        log += info
        log += '\n'

        t2 = datetime.now()
        print('It takes',round((t2-t1).total_seconds()),'seconds to build the channel.')
        print("Start Creating Valley...")
        t1 = datetime.now()
        valley, info = buildValley(paraDict, funDict, channel)
        log += info
        t2 = datetime.now()
        print('It takes',round((t2-t1).total_seconds()),'seconds to build the valley.')
        print('')

        #### Ploting ####
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)

        valleyCol = max(len(channel.levels_x['left']), len(channel.levels_x['right'])) +1
        fig, ax = plt.subplots(1, 1, figsize=[19.2, 14.4], dpi=400)
        ax.plot(channel.x_v*channel.dx, channel.y_center*channel.dx, 'k-', label='CL')
        plotLevels(ax, channel.levels_x, channel.levels_y, channel.dx, 'B', 1)
        plotLevels(ax, valley.levels_x, valley.levels_y, channel.dx, 'V', valleyCol)
        plt.title('SRV Planform')
        plt.xlabel('X (Distance Downstream)')
        plt.ylabel('Y')
        plt.legend()
        #for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
        #             ax.get_xticklabels() + ax.get_yticklabels()):
        #    item.set_fontsize(40)
        plt.savefig(outfolder+'/SRVlevels_xy')

        fig, ax = plt.subplots(1, 1, figsize=[19.2, 14.4], dpi=400)
        ax.plot(channel.x_v*channel.dx, channel.thalweg*channel.dx, 'k-', label='Thalweg')
        plotLevels(ax, channel.levels_x, channel.levels_z, channel.dx, 'B', 1)
        plotLevels(ax, valley.levels_x, valley.levels_z, channel.dx, 'V', valleyCol)
        plt.title('SRV Longitudianl Profile')
        plt.xlabel('X (Distance Downstream)')
        plt.ylabel('Z')
        plt.legend()
        #for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
        #             ax.get_xticklabels() + ax.get_yticklabels()):
        #    item.set_fontsize(40)
        plt.savefig(outfolder+'/SRVlevels_xz')

        fig, ax = plt.subplots(2, 1, sharex=True)
        fig.suptitle('River Centerline Slope & Curvature')
        fig.subplots_adjust(hspace=0)
        plt.subplot(211)
        plt.plot(channel.x_v*channel.dx, channel.getSlope(), 'tab:blue', label='slope')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.legend()
        plt.subplot(212)
        plt.scatter(channel.x_v*channel.dx, channel.getDynamicCurv(), c='tab:blue', s=1, label='dynamic curvature')
        plt.ylabel('Y')
        plt.legend()
        plt.savefig(outfolder+'/SRVcurvature')

        fig = channel.getXShapePlot()
        plt.savefig(outfolder+'/SRVinnerChannelXShape')

        fig = valley.getXShapePlot()
        plt.savefig(outfolder+'/SRVvalleyXShape')
        #### Output files ####

        valley.tocsv(outfolder+"/SRVtopo")
        print("Output files are saved to", os.getcwd()+'/'+outfolder)

        print(log)
        with open(os.getcwd()+'/'+outfolder+'/log.txt', 'w') as f:
            f.write(log)

        print(channel)
        with open(os.getcwd()+'/'+outfolder+'/SRVmetrics.txt', 'w') as f:
            f.write('River Channel Data:\n')
            f.write(channel.__str__())
            f.write('\nValley Data:\n')
            f.write(valley.__str__())

        out = [['X', 'S', 'C']]
        riverSlope = channel.getSlope()
        riverCurvature = channel.getDynamicCurv()
        riverx = channel.x_v*channel.dx
        out += [[riverx[i], riverSlope[i], riverCurvature[i]] for i in range(len(riverx))]
        with open(os.getcwd()+'/'+outfolder+'/SRVcenterline.csv', 'w') as cf:
            wt = csv.writer(cf, lineterminator='\n')
            wt.writerows(out)

        out = [['X', 'Z']]
        xz = valley.tolist_levelxz()
        out += [[xz[0][i]*valley.dx, xz[1][i]*valley.dx] for i in range(len(xz[0]))]
        with open(os.getcwd()+'/'+outfolder+'/SRVlevels_xz.csv', 'w') as cf:
            wt = csv.writer(cf, lineterminator='\n')
            wt.writerows(out)

        out = [['X', 'Y']]
        xz = valley.tolist_levelxy()
        out += [[xz[0][i]*valley.dx, xz[1][i]*valley.dx] for i in range(len(xz[0]))]
        with open(os.getcwd()+'/'+outfolder+'/SRVlevels_xy.csv', 'w') as cf:
            wt = csv.writer(cf, lineterminator='\n')
            wt.writerows(out)

        out = [['X', 'Y', 'ID']]
        xyl = valley.tolabeledlist_levelxy()
        out += [[xyl[0][i]*valley.dx, xyl[1][i]*valley.dx, xyl[2][i]] for i in range(len(xyl[0]))]
        with open(os.getcwd()+'/'+outfolder+'/SRVlevels_xy_labeled.csv', 'w') as cf:
            wt = csv.writer(cf, lineterminator='\n')
            wt.writerows(out)

    except Exception as err:
        print(log)
        traceback.print_exception(*sys.exc_info())
        print(err)
