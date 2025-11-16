import os
# from riverbuilder.core import river
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
import subprocess

np.seterr(divide='ignore', invalid='ignore')

case_names = ['SFE_Leggett_FB_V']

add_thalwegZ = 1
add_width = 1
add_thalwegY = 1
add_valley_wall = 1
add_runway = 1

case_ind = 0
for case_name in case_names:

    gcs_dir = '.\\gcs'
    harmonic_dir = os.path.join('.\\harmonic_functions',
                                case_name)
    case_dir = os.path.join('..\\examples_custom', case_name)
    case_v_dir = os.path.join(case_dir, case_name)

    #############################################################################################
    # Reading gcs excel file - series only
    print('Reading gcs excel file - series only')
    gcs_file = os.path.join(gcs_dir, case_name+'.xlsx')
    gcs = pd.ExcelFile(gcs_file)
    series = pd.read_excel(gcs, 'GVFs')

    topo_file = os.path.join(case_dir, case_name, 'SRVtopo.csv')
    topo_orig = pd.read_csv(topo_file)
    topo = topo_orig

    ind_CL = topo['Label'] == 'CL'
    ind_XS = topo['Label'] == 'XS'
    ind_RB = topo['Label'] == 'RB'
    ind_VL = topo['Label'] == 'VL'

    if add_thalwegZ == 1:
        thalZ_orig = np.max(topo['Y'][topo['Label'] == 'CL'])

        x_CL = topo['X'][ind_CL]
        z_CL_orig = topo['Z'][ind_CL]

        x_fit = series['station']
        z_fit = series['thal_elev'] - (series['thal_elev'][0] - z_CL_orig[0])
        linear_thalZ = interpolate.interp1d(x_fit, z_fit)

        z_CL = linear_thalZ(x_CL)

        x_CL_top = topo.loc[topo['Y'] > 0, 'X'][ind_RB]
        z_CL_top = topo.loc[topo['Y'] > 0, 'Z'][ind_RB]
        linear_CL_top = interpolate.interp1d(x_CL_top, z_CL_top, fill_value='extrapolate') # https://caam37830.github.io/book/04_functions/interpolation.html

        plt.figure()
        plt.plot(x_CL, linear_CL_top(x_CL))
        plt.plot(x_CL, z_CL)
        plt.plot(x_CL, z_CL_orig, 'k')

        sZ = (linear_CL_top(x_CL) - z_CL) / (linear_CL_top(x_CL) - z_CL_orig) ## Using the original slope
        # sZ = (linear_CL_top(x_CL) - z_CL) / (np.max(z_CL_top) - z_CL_orig[0]) ## Adjust river slope = valley slope
        sX = x_CL
        linear_s = interpolate.interp1d(sX, sZ, fill_value='extrapolate')

        plt.figure()
        plt.plot(sX, sZ)

        topo.loc[ind_CL, 'Z'] = z_CL
        x_CL_top_XS = topo['X'][ind_XS]
        z_CL_top_XS = topo['Z'][ind_XS]

        topo.loc[ind_XS, 'Z'] = linear_CL_top(x_CL_top_XS) - \
                                (linear_CL_top(x_CL_top_XS)-z_CL_top_XS) * linear_s(x_CL_top_XS)

        topo_file_thalZ = os.path.join(case_dir, case_name, 'SRVtopo_thalZ1.csv')
        topo.to_csv(topo_file_thalZ)


        plt.figure()
        plt.plot(topo['X'], topo['Y'], '.', color='0.7')
        plt.plot(topo['X'][ind_CL], topo['Y'][ind_CL], 'k.')
        plt.plot(topo['X'][ind_VL], topo['Y'][ind_VL], '.', color='0')
        plt.title(case_name+', Widths added')
        plt.savefig(os.path.join(case_dir, case_name, 'SRVtopo_orig.png'))
        plt.close()

    if add_width == 1:
        width_orig = np.max(topo['Y'][ind_XS])

        ## Left
        x_fit = series['station']
        y_fit = series['l_bankfull'] - series['thal_lat'] * add_thalwegY
        linear_lbf = interpolate.interp1d(x_fit, y_fit)

        dY = linear_lbf(topo['X'][topo['Y'] > 0])
        sY = dY / width_orig
        topo.loc[topo['Y'] > 0, 'Y'] = topo['Y'][topo['Y'] > 0] * sY

        ## Right
        x_fit = series['station']
        y_fit = series['r_bankfull'] - series['thal_lat'] * add_thalwegY
        linear_rbf = interpolate.interp1d(x_fit, y_fit)

        dY = linear_rbf(topo['X'][topo['Y'] < 0])
        sY = -dY / width_orig
        topo.loc[topo['Y'] < 0, 'Y'] = topo['Y'][topo['Y'] < 0] * sY

        ## Valley wall
        if add_valley_wall == 1:
            width_new = np.max(np.abs(topo['Y'][ind_XS]))
            topo.loc[topo['Label'] == 'VL', 'Y'] = np.sign(topo[ind_VL]['Y'])*(width_new + 10)

        plt.figure()
        plt.plot(topo['X'], topo['Y'], '.', color='0.7')
        plt.plot(topo['X'][ind_CL], topo['Y'][ind_CL], 'k.')
        plt.plot(topo['X'][ind_VL], topo['Y'][ind_VL], '.', color='0')
        plt.title(case_name+', Widths added')
        plt.savefig(os.path.join(case_dir, case_name, 'SRVtopo_width.png'))
        plt.close()

    if add_thalwegY == 1:
        x_fit = series['station']
        y_fit = series['thal_lat']
        linear_ThalY = interpolate.interp1d(x_fit, y_fit)

        dY = linear_ThalY(topo['X'])
        topo['Y'] = topo['Y'] + dY

        ## Valley wall
        if add_valley_wall == 1:
            width_new = np.max(np.abs(topo['Y'][topo['Label'] == 'XS']))
            topo.loc[ind_VL, 'Y'] = np.sign(topo[ind_VL]['Y'])*(width_new + 10)

        ## Runway
        if add_runway == 1:
            runway_L = np.max(topo[ind_XS]['X']) / 3
            ind_XS_runway = topo['X'] == topo['X'][0]
            topo.loc[ind_XS_runway, 'X'] = topo[ind_XS_runway]['X'] - runway_L
            ind_XS_runway = topo['X'] == topo['X'][1]
            topo.loc[ind_XS_runway, 'X'] = topo[ind_XS_runway]['X'] - runway_L
            ind_XS_runway = topo['X'] == topo['X'][2]
            topo.loc[ind_XS_runway, 'X'] = topo[ind_XS_runway]['X'] - runway_L

        plt.figure()
        plt.plot(topo['X'], topo['Y'], '.', color='0.7')
        plt.plot(topo['X'][ind_CL], topo['Y'][ind_CL], 'k.')
        plt.plot(topo['X'][ind_VL], topo['Y'][ind_VL], '.', color='0')
        plt.title(case_name+', Thalweg added')
        plt.savefig(os.path.join(case_dir, case_name, 'SRVtopo_thal.png'))
        plt.close()


        plt.figure()
        # plt.plot(topo['X'], topo['Z'], '.')
        plt.plot(topo['X'][ind_XS], topo['Z'][ind_XS], '.', color='0.7')
        plt.plot(topo['X'][ind_CL], topo['Z'][ind_CL], 'k.')
        plt.title(case_name+', Thalweg Z added')
        plt.savefig(os.path.join(case_dir, case_name, 'SRVtopo_thalZ.png'))
        plt.close()


        topo_file_thal = os.path.join(case_dir, case_name, 'SRVtopo_thal.csv')
        topo.to_csv(topo_file_thal)

    case_ind = case_ind + 1

plt.close('all')
