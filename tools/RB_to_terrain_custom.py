#############################################################################################
# RB_to_terrain.py converts the RB output table (.csv) to terrain (.asc)
#   Written by Anzy Lee, Postdoctoral Scholar, Utah State University
#   Date: 10/21/2020
#############################################################################################

import os
import numpy as np
import arcpy
import pandas as pd

#############################################################################################
# Input variables: RB_path, asc_name, RB_unit, asc_unit, cell_size, execute

case_names = ['sfe_tbs_M1_322', 'sfe_tbs_M1_25',
              'sfe_tbs_M1_209', 'sfe_tbs_M1_221',
              'sfe_tbs_M1_24', 'sfe_tbs_M1_81',
              'sfe_tbs_M1_95', 'sfe_tbs_M1_316',
              'sfe_tbs_M1_2248', 'sfe_tbs_M1_82', 'sfe_tbs_M1_4523']

cell_sizes = [1, 1, 0.3, 0.4, 0.7, 1, 1, 1, 0.4, 0.7, 0.2]      # Cell size of ascii terrain in asc_unit
versions = ['vv3', 'vv4']
inner_depth = 'sfe_tbs_M1_bfWSE'


versions = ['vv4']
inner_depth = 'sfe_tbs_M1_halfway'

############################################


case_names = ['SFE_Leggett_yr1_U', 'SFE_Leggett_yr2_U',
              'SFE_Leggett_yr1_V', 'SFE_Leggett_yr2_V',
              'SFE_316_yr1_U', 'SFE_316_yr2_U',
              'SFE_316_yr1_V', 'SFE_316_yr2_V']

inner_depth = 'sfe_tbs_M1_bfWSE'
versions = ['vv4']

cell_sizes = 0.5*np.ones(len(case_names))#[1, 1, 1, 1, 1, 1]

############################################

inner_depth_type = inner_depth.split('_')[-1]

runway = 1
copy_asc_tuflow = 1


RB_unit = 'meter'                # Unit of the river archetype
asc_unit = 'meter'              # Unit of the ascii terrain

case_ind = 0

for case_name in case_names:
    for version in versions:
        case_v_name = case_name + '_' + version

        print('---- '+case_v_name+' ----')

        cell_size = str(np.round(cell_sizes[case_ind],1))
        print('cell size = ' + cell_size)

        execute = np.array([1,          # 1 if you want to execute "Table to point",
                            1,          # 1 if you want to execute "Create TIN",
                            1,          # 1 if you want to execute "TIN to Raster",
                            1])         # 1 if you want to execute "Raster to asc"

        #############################################################################################
        # Workspace setting
        arcpy.env.overwriteOutput = True
        RB_path = os.path.join(os.path.abspath('..'), "samples\\" + inner_depth + "\\"
                               + case_name + "\\" + case_v_name +"\\" + case_v_name)  # path to SRVtopo directory
        arcpy.env.workspace = RB_path
        sr = arcpy.SpatialReference(3857, 115700) #  WGS_1984_web_mercator, WGS 1984
        #sr = arcpy.SpatialReference(4759, 115700) # WGS 1984, WGS 1984
        arcpy.CheckOutExtension("3D")

        if RB_unit == 'foot':
            RB_conv = 3.28084
        else:
            RB_conv = 1
        if asc_unit == 'foot':
            asc_conv = 3.28084
        else:
            asc_conv = 1
        conv_factor = asc_conv/RB_conv

        #############################################################################################
        if execute[0] == 1:
            # 0 Unit conversion
            print('1. Converting units')
            in_Table = arcpy.env.workspace + "\\SRVtopo_thal.csv"
            out_Table = arcpy.env.workspace + "\\SRVtopo_xyz.csv"
            df = pd.read_csv(in_Table)
            offset = 100 # to prevent minus values
            df.X = df.X*conv_factor +offset
            df.Y = df.Y*conv_factor +offset
            df.Z = df.Z*conv_factor +offset
            df.to_csv(out_Table)

            # 1 Table to point
            in_Table = arcpy.env.workspace+"/SRVtopo_xyz.csv"
            output_point = case_v_name+'_xyz_wo_runway.shp'
            x_coords = "X"
            y_coords = "Y"
            z_coords = "Z"

            # Make the XY event layer...
            print("2. Running Table to point")
            arcpy.management.XYTableToPoint(in_Table, output_point,
                                            x_coords, y_coords, z_coords,
                                            sr)
            print("# of points = " + str(arcpy.GetCount_management(output_point)))

            # Points can be adjusted to create a "RUNWAY'
            # print('Points should be adjusted to create a RUNWAY')
            # os.system("pause")

        #############################################################################################
        if execute[1] == 1:
            # 2 Create TIN
            in_point = case_v_name+'_xyz_wo_runway.shp'
            output_TIN = case_v_name+'_wo_runway_TIN'

            print("3. Running Create TIN")
            arcpy.ddd.CreateTin(output_TIN, sr, in_point+" Z masspoints")

        #############################################################################################
        if execute[2] == 1:
            # 3 TIN to Raster
            in_TIN = case_v_name+'_wo_runway_TIN'
            out_tif = case_v_name+'_wo_runway.tif'
            # Set variables for TIN to Raster
            dataType = "FLOAT"  # Default
            method = "LINEAR"  # Default
            sampling = "CELLSIZE " + cell_size
            zfactor = "1"

            print("4. Running TIN Raster")
            arcpy.ddd.TinRaster(in_TIN, out_tif, dataType,
                            method, sampling, zfactor)
        #############################################################################################

        if runway == 1:
            in_raster = case_v_name + '_wo_runway.tif'
            out_point_features = case_v_name + '_xyz.shp'
            arcpy.conversion.RasterToPoint(in_raster, out_point_features)

            # input("5. Press enter after making a runway using "+out_point_features)

            print("---- converting point to TIN ----")
            in_point = case_v_name + '_xyz.shp'
            out_TIN = case_v_name + '_TIN'
            arcpy.ddd.CreateTin(out_TIN, sr, in_point + " grid_code masspoints")

            print("---- converting TIN to raster ----")
            in_TIN = case_v_name + '_TIN'
            out_tif = case_v_name + '.tif'
        else:
            print("---- converting TIN to raster ----")
            in_TIN = case_v_name + '_TIN_wo_runway'
            out_tif = case_v_name + '.tif'

        # Set variables for TIN to Raster
        dataType = "FLOAT"  # Default
        method = "LINEAR"  # Default
        sampling = "CELLSIZE " + cell_size
        zfactor = "1"


        arcpy.ddd.TinRaster(in_TIN, out_tif, dataType,
                            method, sampling, zfactor)

        if execute[3] == 1:
            # 4 Raster to ascii
            in_tif = case_v_name+'.tif'
            out_ascii = os.path.join(RB_path, case_v_name + '.asc')

            print("6. Running Raster To ASCII")
            arcpy.RasterToASCII_conversion(in_tif, out_ascii)

        if copy_asc_tuflow == 1:

            terrain_dir = 'E:\\tuflow-RB-tbs\\'+inner_depth_type+'\\terrains\\' + case_v_name
            if not os.path.exists(terrain_dir):
                os.mkdir(terrain_dir)

            in_raster = out_ascii
            out_rasterdataset = os.path.join(terrain_dir, case_v_name + '.asc')
            arcpy.management.Copy(in_raster, out_rasterdataset)
        #############################################################################################

    case_ind = case_ind + 1