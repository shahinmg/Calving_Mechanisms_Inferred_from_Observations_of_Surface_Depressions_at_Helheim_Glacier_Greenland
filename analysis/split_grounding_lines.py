#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 11:57:28 2025

@author: laserglaciers
"""

import geopandas as gpd
import os

path_2019 = '/media/laserglaciers/upernavik/flexure_manuscipt/zenodo_submission/data/ground_lines_and_grounding_zones/grounding_lines_merged/2018_grounding_lines_merge.gpkg'


df2019 = gpd.read_file(path_2019)


date_1_set = set(df2019['date_1'])

date_row_list = []
multi_dates = []
date_len_sum = 0
for date in date_1_set:
    
    
    singedate_msk = df2019['date_1'] == date
    single_date = df2019[singedate_msk]
    
    if len(single_date) == 1:
        date_row_list.append(single_date)
    
    elif len(single_date) > 1:
        
        date_len_sum += len(single_date)
        multi_dates.append(date)
        
        multi_line = single_date.dissolve()
        date_row_list.append(multi_line)

df_fix = gpd.pd.concat(date_row_list)

multi_line_df = df2019[df2019['date_1'].isin(multi_dates)]

op = '/media/laserglaciers/upernavik/flexure_manuscipt/zenodo_submission/data/'

df_fix.to_file(f'{op}2018_grounding_lines_merged_v2.gpkg')
