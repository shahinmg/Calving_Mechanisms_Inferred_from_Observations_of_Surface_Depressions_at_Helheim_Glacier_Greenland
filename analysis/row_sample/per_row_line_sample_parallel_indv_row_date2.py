#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 12:35:45 2023

@author: laserglaciers
"""

import os
import pandas as pd
import geopandas as gpd
import rasterio
import numpy as np
from line_sample_tool import *
from subProcessUtils import date_sort,date_sort_two_dates
from multiprocessing import Pool, cpu_count


row_line_path = '/media/laserglaciers/QAANAAQ/atlas_south_data_full_catalog/individual_rows_extended/'
data_path = '/media/laserglaciers/upernavik/atlas_south_copc/velocity_standard_3413/smoothed/distance8/all_bands/'
op = '/media/laserglaciers/upernavik/atlas_south_copc/pkls/abs_velocity_standard_d8/abs_velocity_standard_d8_date2/'


row_list = sorted([pkl for pkl in os.listdir(row_line_path) if pkl.endswith('gpkg')],
                  key=lambda x:int(x.split('_')[-1].split('.')[0]))

row_list = [row for row in os.listdir(row_line_path) if row.endswith('gpkg')]
dem_tups, dem_files, dem_times2 = date_sort_two_dates(data_path, 'tif')

date1 = np.datetime64('2015')
date2 = np.datetime64('2024')
date_tups_filtered = [(f'{data_path}{file}',dates2) for file,dates2 in dem_tups if dates2>=date1 and dates2<=date2]


arg_list = []
os.chdir(data_path)
for row, line in enumerate(row_list):
    
    line_path = f'{row_line_path}{line}'
    op_final = f'{op}row_{row}_abs_velocity_standard_d8_date2'
    arg_list.append((line_path,date_tups_filtered[:],op_final))
    # sample_points_linestring_from_rasters_v2(line,dem_tups,op_final,dem=False,gpd_read_in=False) # DEM = False bc these are single band dems
    
# MAKE SURE TO CHANGE WHICH BAND TO SAMPLE FOR DEM OR VELOCITY
if __name__ == "__main__":
    with Pool(processes=cpu_count()-2) as pool:
          results = pool.starmap(sample_points_linestring_from_rasters_v3_par, arg_list)
