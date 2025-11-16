
import os
# from riverbuilder.core import river
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fRB import *
from signal_to_riverbuilder_detrend import *
import subprocess

case_names = ['SFE_Leggett_FB_V']

avg_bfw, channel_slope = [], []
exec_harmonic = 1

for case_name in case_names:

    gcs_dir = '.\\gcs'
    harmonic_dir = os.path.join('.\\harmonic_functions',
                                case_name)
    case_dir = os.path.join('..\\examples_custom', case_name)

    #############################################################################################
    # Reading gcs excel file - series only
    print('Reading gcs excel file - series only')
    gcs_file = os.path.join(gcs_dir, case_name+'.xlsx')
    gcs = pd.ExcelFile(gcs_file)
    series = pd.read_excel(gcs, 'GVFs')
    Z_d, valley_slope = get_detrended(series)

    variables = ['Z_d', 'W_bf_L', 'W_bf_R', 'thal_y']

    dataset = np.transpose([series.station, Z_d,
                            series.l_bankfull - series.thal_lat,
                            series.thal_lat - series.r_bankfull,
                            series.thal_lat])

    # Writing an excel file for harmonic analysis
    df = pd.DataFrame(dataset, index=None, columns=['Station'] + variables)
    blankIndex = [''] * len(df)
    df.index = blankIndex
    harmonic_file = os.path.join(harmonic_dir, case_name + '.csv')

    if not os.path.exists(harmonic_dir):
        os.mkdir(harmonic_dir)
    df.to_csv(harmonic_file)
    #############################################################################################

    if exec_harmonic == 1:
        in_csv = harmonic_file
        lat_offset, slope, y_interc = river_builder_harmonics(in_csv, 'Station', units='m', fields=variables, field_names=[],
                                             r_2=0.99, n=20, methods='by_fft')
        # lat_offset = minimum lateral offset of variables in order of "variables'
        # matplotlib.pyplot.close('all')

    geo_params = ini_geo_params()

    # min_W_base = lat_offset[1] * 2
    # min_W_bf = lat_offset[1] * 2
    min_W_bf_L = lat_offset[1]
    min_W_bf_R = lat_offset[2]
    geo_params = write_geo_params(geo_params, 'min_W_bf_L', min_W_bf_L)
    geo_params = write_geo_params(geo_params, 'min_W_bf_R', min_W_bf_R)


    #############################################################################################
    # Making Harmonic function txt input files for RB
    if not os.path.exists(case_dir):
        os.mkdir(case_dir)

    ind_v = 0
    for variable in variables:
        print(variable)
        in_csv = os.path.join(harmonic_dir, variable + '_harmonics_by_fft.csv')
        out_txt = os.path.join(case_dir, variable + '.txt')
        stat = "python harmonicParser.py " + in_csv + " " + out_txt
        print(stat)
        p = subprocess.run(stat)


    #############################################################################################
    # Calculating the metrics in gcs excel file
    geo_params = get_geo_params(gcs_file, geo_params, case_name, 'vv4')

    #############################################################################################
    # Writing RB input txt file
    f = open(os.path.join(case_dir, case_name + '.txt'), 'w+')
    RB_params = ini_RB_params()

    RB_params = version_params('vv4', RB_params, geo_params)

    if case_name.split('_')[-1] == 'V':
        RB_params['Xshape'][0] = 'V'

    disclaimer(f)
    domain_params(f, RB_params)
    channel_params(f, RB_params)
    crossX(f, RB_params)
    userfuncs(f, RB_params)
    channelbreakparams(f, RB_params)
    # bedroughness(f, RB_params) ## if reach length is too big, does not work sometimes
    f.close()

    avg_bfw = np.append(avg_bfw, np.average(series.l_bankfull*2))
    channel_slope = np.append(channel_slope, RB_params['valley_slope'][0])

channel_att = pd.DataFrame([case_names, avg_bfw, channel_slope]).transpose()
channel_att.columns = ['case_names', 'avg_bfw', 'slope']
channel_att['avg_bfw'] = channel_att['avg_bfw'].astype('float')
channel_att['slope'] = channel_att['slope'].astype('float')