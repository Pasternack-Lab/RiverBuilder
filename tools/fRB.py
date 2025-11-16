import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def get_detrended(series):
    # Valley slope
    # - Linear regression
    Station_d1 = np.vstack(series.station)  # adding a dimension
    X = np.hstack((np.ones((np.size(Station_d1), 1)), Station_d1))
    Y = np.vstack(series.thal_elev)
    A = np.linalg.inv(X.transpose().dot(X)).dot(X.transpose()).dot(Y)
    valley_slope = A[1][0]
    y_interc = A[0][0]
    Z_d = series.thal_elev - (valley_slope * series.station + y_interc)
    return Z_d, valley_slope

def ini_geo_params():
    geo_params = np.array([(0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 'SU')],
                          dtype=[('min_W_base', 'f4'), ('min_h_base', 'f4'), ('avg_h_bf', 'f4'),
                                 ('min_W_bf', 'f4'), ('min_W_bf_L', 'f4'), ('min_W_bf_R', 'f4'), ('min_h_bf', 'f4'),
                                 ('avg_W_bf', 'f4'), ('domain_length', 'f4'), ('amp', 'f4'),
                                 ('freq', 'f4'), ('phase', 'f4'),
                                 ('min_W_TIN', 'f4'), ('avg_h_TIN', 'f4'),
                                 ('valley_slope', 'f4'), ('PBR', 'f4'), ('Xshape', 'U4')]
                          )
    return geo_params

def ini_RB_params():
    RB_params = np.array([(0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 'SU', 'vv0', 0, 0.5)],
                         dtype=[('min_inner_lat', 'f4'), ('min_inner_depth', 'f4'), ('domain_length', 'f4'),
                                ('valley_slope', 'f4'),
                                ('amp', 'f4'), ('freq', 'f4'), ('phase', 'f4'),
                                ('Zd_m', 'U20'), ('L_inner', 'U20'), ('R_inner', 'U20'),
                                ('L1_outer_lat_min', 'f4'), ('L1_outer_h', 'f4'), ('L1_outer_func', 'U20'),
                                ('R1_outer_lat_min', 'f4'), ('R1_outer_h', 'f4'), ('R1_outer_func', 'U20'),
                                ('L2_outer_lat_min', 'f4'), ('L2_outer_h', 'f4'), ('R2_outer_lat_min', 'f4'),
                                ('R2_outer_h', 'f4'),
                                ('PBR', 'f4'), ('Xshape', 'U10'), ('version', 'U20'), ('valley_height_offset', 'f4'),
                                ('X_resolution', 'f4')])
    return RB_params

def write_geo_params(geo_params, field, value):
    geo_params[field] = value
    return geo_params

def get_geo_params(gcs_file, geo_params, case_name, version):
    gcs = pd.ExcelFile(gcs_file)

    # if case_name.split('_')[1] == 'tbs' or case_name[0] == 'M':
    #     series = pd.read_excel(gcs, 'GVFs')
    # else:
    #     series = pd.read_excel(gcs, 'series')

    series = pd.read_excel(gcs, 'GVFs')
    # base_bf_metrics = pd.read_excel(gcs, 'base_bf')  # having both base & bf width var.
    # bf_metrics = pd.read_excel(gcs, 'bf')  # having only bf width var.

    #############################################################################################
    # Calculating metrics

    #geo_params['min_h_base'] = np.min(series.WSE_base - series.Z)
    geo_params['avg_h_bf'] = np.average(series.bkf_elev - series.thal_elev)
    geo_params['min_h_bf'] = np.min(series.bkf_elev - series.thal_elev)

    if version in ['vv1', 'vv3', 'vv4']:
        print(version)
        # geo_params['domain_length'] = base_bf_metrics.Value[0]
        geo_params['domain_length'] = max(series.station)
        geo_params['amp'] = 0 #base_bf_metrics.Value[3]
        geo_params['freq'] = 0 #base_bf_metrics.Value[4]
        geo_params['phase'] = 0 #base_bf_metrics.Value[5]
        geo_params['min_W_TIN'] = 0 #base_bf_metrics.Value[15]
        geo_params['avg_h_TIN'] = 0 #base_bf_metrics.Value[16]
        geo_params['avg_W_bf'] = np.average(series.l_bankfull*2)
        geo_params['PBR'] = 0 ## vv1 and vv3
        geo_params['Xshape'] = 'AU'
        #geo_params['Xshape'] = base_bf_metrics.Value[25]

        if case_name.split('_')[2] == 'M1':
            geo_params['domain_length'] = max(series.station)
            geo_params['amp'] = 0
            geo_params['freq'] = 0
            geo_params['phase'] = 0
            geo_params['min_W_TIN'] = max(max(series.l_bankfull), max(series.r_bankfull)) + 5

            if version == 'vv4':
                geo_params['PBR'] = 0#base_bf_metrics.Value[24] # vv4

        if version in ['vv3', 'vv4']:
            geo_params['min_W_bf'] = np.min(series.l_bankfull-series.r_bankfull)
    else:
        geo_params['domain_length'] = max(series.Station)
        if case_name.split('_')[1] == 'tbs':
            geo_params['domain_length'] = max(series.Station)
        geo_params['amp'] = base_bf_metrics.Value[3]
        geo_params['freq'] = base_bf_metrics.Value[4]
        geo_params['phase'] = base_bf_metrics.Value[5]
        geo_params['min_W_TIN'] = base_bf_metrics.Value[15]
        geo_params['avg_h_TIN'] = base_bf_metrics.Value[16]
        geo_params['avg_W_bf'] = np.average(series.W_bf)
        geo_params['PBR'] = base_bf_metrics.Value[24]
        geo_params['Xshape'] = base_bf_metrics.Value[25]

    Z_d, valley_slope = get_detrended(series)

    geo_params['valley_slope'] = valley_slope
    return geo_params

def version_params(version, RB_params, geo_params):
    Zd_m = 'Z_d.txt'
    RB_params['domain_length'] = geo_params['domain_length']
    RB_params['valley_slope'] = geo_params['valley_slope']
    RB_params['amp'] = geo_params['amp']
    RB_params['freq'] = geo_params['freq']
    RB_params['phase'] = geo_params['phase']
    RB_params['Zd_m'] = Zd_m
    RB_params['PBR'] = geo_params['PBR']
    RB_params['Xshape'] = geo_params['Xshape'][0]
    RB_params['version'] = version

    if version == 'vv0' or version == 'r0':
        L_inner = 'W_base_half.txt'
        R_inner = 'W_base_half.txt'
        L1_outer_func = 'W_bf_half.txt'
        R1_outer_func = 'W_bf_half.txt'

        RB_params['min_inner_lat'] = geo_params['min_W_base']
        RB_params['min_inner_depth'] = geo_params['min_h_base']
        RB_params['L_inner'] = L_inner
        RB_params['R_inner'] = R_inner
        RB_params['L1_outer_lat_min'] = geo_params['min_W_bf']/2
        RB_params['L1_outer_h'] = geo_params['avg_h_bf']
        RB_params['L1_outer_func'] = L1_outer_func
        RB_params['R1_outer_lat_min'] = geo_params['min_W_bf']/2
        RB_params['R1_outer_h'] = geo_params['avg_h_bf']
        RB_params['R1_outer_func'] = R1_outer_func
        RB_params['L2_outer_lat_min'] = geo_params['min_W_TIN']
        RB_params['L2_outer_h'] = geo_params['avg_h_TIN']
        RB_params['R2_outer_lat_min'] = geo_params['min_W_TIN']
        RB_params['R2_outer_h'] = geo_params['avg_h_TIN']

        if version == 'r0':
            RB_params['PBR'] = 0

    else:

        L_inner = 'W_bf_half.txt'
        R_inner = 'W_bf_half.txt'

        RB_params['min_inner_lat'] = geo_params['min_W_bf']
        RB_params['min_inner_depth'] = geo_params['min_h_bf']
        RB_params['L_inner'] = L_inner
        RB_params['R_inner'] = R_inner
        RB_params['L1_outer_lat_min'] = geo_params['min_W_TIN']
        RB_params['L1_outer_h'] = geo_params['avg_h_TIN']
        RB_params['R1_outer_lat_min'] = geo_params['min_W_TIN']
        RB_params['R1_outer_h'] = geo_params['avg_h_TIN']

        if version == 'c0':
            RB_params['PBR'] = 0
        elif version == 'c1':
            RB_params['L_inner'] = '0'
            RB_params['R_inner'] = '0'
            RB_params['min_inner_lat']=geo_params['avg_W_bf']
        elif version == 'c2':
            RB_params['Zd_m'] = '0'
            RB_params['min_inner_depth'] = geo_params['avg_h_bf']
        elif version == 's0':
            RB_params['L_inner'] = '0'
            RB_params['R_inner'] = '0'
            RB_params['min_inner_lat'] = geo_params['avg_W_bf']
            RB_params['min_inner_depth'] = geo_params['avg_h_bf']
            RB_params['Zd_m'] = '0'
        elif version == 's1':
            RB_params['PBR'] = 0
            RB_params['Zd_m'] = '0'
            RB_params['min_inner_depth'] = geo_params['avg_h_bf']
        elif version == 's2':
            RB_params['L_inner'] = '0'
            RB_params['R_inner'] = '0'
            RB_params['min_inner_lat'] = geo_params['avg_W_bf']
            RB_params['PBR'] = 0
        elif version == 'vv1':
            RB_params = RB_params # all the parameters but baseflow width variation
        elif version in ['vv3', 'vv4']:
            RB_params['L_inner'] = 'W_bf_L.txt'
            RB_params['R_inner'] = 'W_bf_R.txt' # asymmetric
            # RB_params['R_inner'] = 'W_bf_L.txt' #symmetric

            RB_params['min_inner_lat'] =geo_params['min_W_bf']
        elif version == 'r1':
            RB_params['min_inner_depth'] = 0.001
            RB_params['Zd_m'] = '0'
            RB_params['valley_slope'] = 0
        elif version == 'n0':
            RB_params['L_inner'] = '0'
            RB_params['R_inner'] = '0'
            RB_params['min_inner_lat'] = geo_params['avg_W_bf']
            RB_params['min_inner_depth'] = geo_params['avg_h_bf']
            RB_params['Zd_m'] = '0'
            RB_params['PBR'] = 0
        else:
            print('ERROR: You can only choose the following versions: vv0, vv1, vv3, r0, r1, c0, c1, c2')

    return RB_params

def disclaimer(f):
    f.write("\n#DISCLAIMER" +
            "\n\n#No warranty is expressed or implied regarding the usefulness or completeness of "
            "the information provided by River Builder and its documentation. "
            "References to commercial products do not imply endorsement by the Authors of River Builder. "
            "The concepts, materials, and methods used in the algorithms and described in the manual "
            "are for informational purposes only. The Authors have made substantial effort to ensure "
            "the accuracy of the algorithms and the manual, but science is uncertain and the Authors "
            "nor their employers or funding sponsors shall not be held liable for calculations and/or "
            "decisions made on the basis of application of River Builder. The information is provided \"as is\" "
            "and anyone who chooses to use the information is responsible for her or his own choices as to "
            "what to do with the data and the individual is responsible for the results the follow from their "
            "decisions." +
            "\n\n# This input is intended to be used for riverbuilder1.0.0." +
            "\n\n###############" + "\n\n#### NOTES ####" + "\n\n###############" +
            "\n\n# - All dimensional numbers are in units of meters." +
            "\n\n# - User-defined functions may be used for sub-reach variability parameters only." +
            "\n\n# - Everyline starting with '#' will be ignored. If one wants to use example inputs, s/he needs to" +
            "\n\n#   delete the '#' at the start of paramenter lines." +
            "\n\n# - Bankfull depth can either be (A) user-defined or (B) calculated from the Critical Shields Stress "
            "and Median Sediment Size." +
            "\n\n# - Centerline Curvature can either be (A) user-defined or (B) calculate from centerline slope." +
            "\n\n# - Calculations of banks of channel are based on channel centerline; calculations of levels of valley "
            "are based on valley centerline." +
            "\n\n####################" + "\n\n#### GUIDELINES ####" + "\n\n####################" +
            "\n\n#   To ensure the program functions correctly, please abide by the following:" +
            "\n\n#   - Put each parameter or function in a seperate line." +
            "\n\n#   - When providing an input involving pi, do so in the form of a*pi, where a is a constant. "
            "(EX: 2*pi, pi, pi/6, 3+pi)"
            )

def domain_params(f, RB_params):
    f.write("\n\n####################################" + "\n\n#### DOMAIN PARAMETERS (METERS) ####" + "\n\n####################################" +
            "\n\n# These parameters initial the environment of the river; the follow will be an example input." +
            "\n\n# The \"Smooth\" variable prevents sharp turning in centerline by averaging adjacent values." +
            "\n#    Smooth=1 means taking average from 1 left point and 1 right point, so it will be 3 points in" +
            "\n#    total (including itself). Smooth=0 means no smooth at all. Whenever piece-wise mask function" +
            "\n#    is applied on centerline, at least Smooth=1 is recommended." +
            "\n\nDatum=500" +
            "\nLength=%9.3f" % RB_params['domain_length'] +
            "\nValley Slope (Sv)=%.4f" % -np.round(RB_params['valley_slope'],4) +
            "\nSmooth=1" +
            "\nX Resolution=%.1f" % RB_params['X_resolution']
            )

def channel_params(f, RB_params):
    f.write("\n\n####################################" + "\n\n#### CHANNEL PARAMETERS (METERS) ####" + "\n\n####################################" +
            "\n\n# These parameters controls the inner channel of the river; the follow will be an example input." +
            "\n# Note that if \"Inner Channel Depth Minimum\" is set to 0, then the depth minimum will be " +
            "\n# calculated based on median sediment size and critical shields stress." +
            "\n\nInner Channel Lateral Offset Minimum=%.4f" % RB_params['min_inner_lat'] +
            "\nInner Channel Depth Minimum=%.4f" % RB_params['min_inner_depth'] +
            "\nMedian Sediment Size (D50)=0" +
            "\nCritical Shields Stress (t*50)=0"
            )

def crossX(f, RB_params):
    if RB_params['Xshape'][0] == 'AU':
        f.write("\n\n###############################" + "\n\n#### CROSS SECTIONAL SHAPE ####" + "\n\n###############################" +
                "\n\n# These parameters handles the cross-sectional shape of inner channel. " +
                "\n\nChannel XS Points=21" +
                "\nCross-Sectional Shape=" + RB_params['Xshape'][0]
                )
    elif RB_params['Xshape'][0] == 'V':
        f.write("\n\n###############################" + "\n\n#### CROSS SECTIONAL SHAPE ####" + "\n\n###############################" +
                "\n\n# These parameters handles the cross-sectional shape of inner channel. " +
                "\n\nChannel XS Points=21" +
                "\nCross-Sectional Shape=EN" +
                "\nTZ(n)=1"
                )

def userfuncs(f, RB_params):
    f.write("\n\n################################" + "\n\n#### USER-DEFINED FUNCTIONS ####" + "\n\n################################" +
            "\n\nMASK0=(ALL)" +
            "\nSIN1=(%.4f, %.4f, %.4f, MASK0)" % (RB_params['amp'], RB_params['freq'], RB_params['phase'])
            )

def channelbreakparams(f, RB_params):
    if RB_params['version'][0] == 'vv0' or RB_params['version'][0] == 'r0':
        f.write("\n\n############################################" + "\n\n#### CHANNEL BREAKLINE INPUT PARAMETERS ####" + "\n\n############################################" +
                "\n\nMeandering Centerline Function=SIN1" +
                "\nThalweg Elevation Function=" + RB_params['Zd_m'][0] +
                "\n\n# Inner Bank Properties" +
                "\nLeft Inner Bank Function=" + RB_params['L_inner'][0] +
                "\nRight Inner Bank Function=" + RB_params['R_inner'][0] +
                "\n\n# Left Outer Banks Properties" +
                "\nL1 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['L1_outer_lat_min'] +
                "\nL1 Outer Bank Height Offset=%.4f" % RB_params['L1_outer_h'] +
                "\nL1 Outer Bank Function=" + RB_params['L1_outer_func'][0] +
                "\n\n# Right Outer Banks Properties" +
                "\nR1 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['R1_outer_lat_min'] +
                "\nR1 Outer Bank Height Offset=%.4f" % RB_params['R1_outer_h'] +
                "\nR1 Outer Bank Function=" + RB_params['R1_outer_func'][0] +
                "\n\n# Second outer banks" +
                "\n\n# Left Outer Banks Properties" +
                "\nL2 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['L2_outer_lat_min'] +
                "\nL2 Outer Bank Height Offset=%.4f" % RB_params['L2_outer_h'] +
                "\n\n# Right Outer Banks Properties" +
                "\nR2 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['R2_outer_lat_min'] +
                "\nR2 Outer Bank Height Offset=%.4f" % RB_params['R2_outer_h'] +
                "\n\n# Left Valley Boundary Properties" +
                "\nLeft Valley Boundary Lateral Offset Minimum=10" +
                "\nLeft Valley Boundary Height Offset=%.4f" % RB_params['valley_height_offset'] +
                "\n\n# Right Valley Boundary Properties" +
                "\nRight Valley Boundary Lateral Offset Minimum=10" +
                "\nRight Valley Boundary Height Offset=%.4f" % RB_params['valley_height_offset']
                )
    elif RB_params['version'][0] in ['vv3', 'vv4']:
        f.write("\n\n############################################" + "\n\n#### CHANNEL BREAKLINE INPUT PARAMETERS ####" + "\n\n############################################" +
                # "\n\nMeandering Centerline Function=thal_y.txt" +
                # "\nThalweg Elevation Function=" + RB_params['Zd_m'][0] +
                # "\n\n# Inner Bank Properties" +
                # "\nLeft Inner Bank Function=" + RB_params['L_inner'][0] +
                # "\nRight Inner Bank Function=" + RB_params['R_inner'][0] +
                # "\n\n# Left Outer Banks Properties" +
                # "\nL1 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['L1_outer_lat_min'] +
                # "\nL1 Outer Bank Height Offset=%.4f" % RB_params['L1_outer_h'] +
                # "\nL1 Outer Bank Function=" + RB_params['L1_outer_func'][0] +
                # "\n\n# Right Outer Banks Properties" +
                # "\nR1 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['R1_outer_lat_min'] +
                # "\nR1 Outer Bank Height Offset=%.4f" % RB_params['R1_outer_h'] +
                # "\nR1 Outer Bank Function=" + RB_params['R1_outer_func'][0] +
                "\n\n# Left Valley Boundary Properties" +
                "\nLeft Valley Boundary Lateral Offset Minimum=10" +
                "\nLeft Valley Boundary Height Offset=%.4f" % RB_params['valley_height_offset'] +
                "\n\n# Right Valley Boundary Properties" +
                "\nRight Valley Boundary Lateral Offset Minimum=10" +
                "\nRight Valley Boundary Height Offset=%.4f" % RB_params['valley_height_offset']
                )
    else:
        f.write("\n\n############################################" + "\n\n#### CHANNEL BREAKLINE INPUT PARAMETERS ####" + "\n\n############################################" +
                "\n\nMeandering Centerline Function=SIN1" +
                "\nThalweg Elevation Function=" + RB_params['Zd_m'][0] +
                "\n\n# Inner Bank Properties" +
                "\nLeft Inner Bank Function=" + RB_params['L_inner'][0] +
                "\nRight Inner Bank Function=" + RB_params['R_inner'][0] +
                "\n\n# Left Outer Banks Properties" +
                "\nL1 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['L1_outer_lat_min'] +
                "\nL1 Outer Bank Height Offset=%.4f" % RB_params['L1_outer_h'] +
                "\n\n# Right Outer Banks Properties" +
                "\nR1 Outer Bank Lateral Offset Minimum=%.4f" % RB_params['R1_outer_lat_min'] +
                "\nR1 Outer Bank Height Offset=%.4f" % RB_params['R1_outer_h'] +
                "\n\n# Left Valley Boundary Properties" +
                "\nLeft Valley Boundary Lateral Offset Minimum=10" +
                "\nLeft Valley Boundary Height Offset=5" +
                "\n\n# Right Valley Boundary Properties" +
                "\nRight Valley Boundary Lateral Offset Minimum=10" +
                "\nRight Valley Boundary Height Offset=5"
                )

def bedroughness(f, RB_params):

    f.write("\n\n###################################" + "\n\n########## Bed Roughness ##########" + "\n\n###################################" +
            "\n\n# PBR - Perlin bed roughness. Has only one height attribute. Therefore the noise ranges " +
            "\n#        from (-height, height)."
            "\n\nPBR=%.1f" % RB_params['PBR']
            )

