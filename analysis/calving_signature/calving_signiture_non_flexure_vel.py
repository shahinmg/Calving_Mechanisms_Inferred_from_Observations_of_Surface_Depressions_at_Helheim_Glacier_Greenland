#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:27:20 2024

@author: laserglaciers
"""

# =============================================================================
# Need to get the center of the calving block 
# Find the closest profile line and sample that
# I should do this for both flexure and non-flexure calving events
# I should maybe try to sample the full with of the calving bit?
# Though that is probably too much and not enough will happen.
# =============================================================================

import os
import pandas as pd
import geopandas as gpd
import numpy as np
from subProcessUtils import date_sort
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm, colors
from pathlib import Path
from shapely.geometry import Point
import numpy as np
from matplotlib import colors, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

vel_d8_pkl_path = '/media/laserglaciers/upernavik/atlas_south_copc/pkls/abs_velocity_d16_24hr/'
exx_st_d8_pkl_path = '/media/laserglaciers/upernavik/atlas_south_copc/pkls/strain_rates_d8/'
calving_blocks_path = '/media/laserglaciers/QAANAAQ/flexure_vs_non-flexure_calving/merged_calving_blocks/combined_calving_blocks_volume.gpkg'
tides_path = '/media/laserglaciers/upernavik/atlas_south_copc/tide_pkl/tmd_modeled_tide_all_components-2015-2023.pkl'
row_paths = '/media/laserglaciers/upernavik/atlas_sample_rows/merged_rows.gpkg'


vel_d8_pkl_list = sorted([pkl for pkl in os.listdir(vel_d8_pkl_path) if pkl.endswith('pkl')],
                  key=lambda x:int(x.split('_')[1]))


exx_st_d8_pkl_list = sorted([pkl for pkl in os.listdir(exx_st_d8_pkl_path) if pkl.endswith('pkl')],
                  key=lambda x:int(x.split('_')[1]))

# profile_rows = sorted([gpkg for gpkg in os.listdir(row_paths) if gpkg.endswith('gpkg')],
#                   key=lambda x:int(x.split('_')[1][:-5]))
profile_rows = gpd.read_file(row_paths)

flexure_data_path = '/media/laserglaciers/upernavik/atlas_south_copc/calving_catalog_2015-2023/flexure_vs_non-flexure_calving_catalog.csv'
flexure_data = pd.read_csv(flexure_data_path ,parse_dates=(['first_scan', 'last_scan']))

# get dates based on calving style
flexure = flexure_data[flexure_data['type'] == 'flexure']
non_flexure = flexure_data[flexure_data['type'] == 'non-flexure']
flexure_dates = flexure['first_scan']
non_flexure_dates = non_flexure['first_scan']


dem_max_path = '/media/laserglaciers/upernavik/atlas_south_copc/pkls/dem_profiles/max_geoid_rm/'
term_path = '/media/laserglaciers/upernavik/atlas_terminus_finder/terminus_position_elevation_dfs/'
z_displacement_path = '/media/laserglaciers/upernavik/atlas_south_copc/pkls/z_displacement/z_displacement_date2_d3/'

dem_max_pkl_list = sorted([pkl for pkl in os.listdir(dem_max_path) if pkl.endswith('pkl')],
                  key=lambda x:int(x.split('_')[1]))

term_gdf_list = sorted([pkl for pkl in os.listdir(term_path) if pkl.endswith('gpkg')],
                  key=lambda x:int(x.split('_')[1].split('.')[0]))

z_disp_pkl_list = sorted([pkl for pkl in os.listdir(z_displacement_path) if pkl.endswith('pkl')],
                  key=lambda x:int(x.split('_')[1]))

tide = pd.read_pickle(tides_path)
calving_blocks = gpd.read_file(calving_blocks_path)



def multi_row_avergage(pkl_num_list,data_list,term_list,first_scan_date,pkl_path):
    
    avg_list = []
    avg_term_pos_list = []
    avg_elevation_list = []
    for pkl_num in pkl_num_list:
        
        # os.chdir(dem_surface_path)
        # term_df = pd.read_pickle(term_pkl_list[pkl_num])
        # coords = term_df['coords']

        # term_df = term_df.iloc[:,:-2]
        # term_df.columns = pd.Datfrom matplotlib import cmetimeIndex(term_df.columns)
        # term_df = term_df.where(term_df != 0.0, np.nan)
        
        os.chdir(pkl_path)
        pkl_df = pd.read_pickle(data_list[pkl_num])
        pkl_df = pkl_df.iloc[:,:-2]
        pkl_df.columns = pd.DatetimeIndex(pkl_df.columns)
        
    return 

# I need to find which profile line is closest to the center of the calving block
# First get the centroid of the calving block
center_profile_calving_block = {}
for index, row in calving_blocks.iterrows():
    if row['geometry'] != None:
        centroid = row['geometry'].centroid # center of mass assuming constant density
        
    # # Need to find which line is closest to the centroid
    distance_dict = {}
    for p_idx, profile_row in profile_rows.iterrows():
        geom = profile_row['geometry']
        distance_dict[profile_row['layer']] = centroid.distance(geom)
    
    # create pandas series from dictionary and get idxmin
    distance_series = pd.Series(distance_dict)
    center_profile_calving_block[index] = distance_series.idxmin()
    center_profile_calving_block_series = pd.Series(center_profile_calving_block)
    

calving_blocks['center_row'] = center_profile_calving_block_series
flexure_blocks = calving_blocks[calving_blocks['type']=='flexure']
non_flexure_blocks = calving_blocks[calving_blocks['type']=='non-flexure']

fontsize = 20
x_lim = (301548.2189, 314088.2189)
title_list = ['$\dot{\epsilon}_{xx}$ 6hr d8',
              # '$\dot{\epsilon}_{xx}$ 6hr d3',
               # '$\dot{\epsilon}_{xx}$ 24hr d8',
              # '$\dot{\epsilon}_{xx}$ 24hr d3'
              ]
#%% 
# plot flexure calving 
# for date1, date2 in zip(non_flexure['first_scan'], non_flexure['last_scan']):
op = '/media/laserglaciers/upernavik/atlas_south_copc/figs/flexure_vs_non_flexure/behind_calving_block/non-flexure/individual_calving/'
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,8))
# average_vel_arr = np.empty((56,17))
# average_vel_arr = average_vel_arr[:] * np.nan
# average_day_arr = average_vel_arr.copy()

df_store = []

i = 0
for block_idx, block_row in non_flexure_blocks.iterrows():
    

    date1 = block_row['first_date']
    date2 = block_row['second_date']
    profile_row = block_row['center_row']
    row_number = int(profile_row.split('_')[-1])
    
    vel_pkl_full = pd.read_pickle(f'{vel_d8_pkl_path}{vel_d8_pkl_list[row_number]}')
    exx_pkl = pd.read_pickle(f'{exx_st_d8_pkl_path}{exx_st_d8_pkl_list[row_number]}')
    term_gdf = gpd.read_file(f'{term_path}{term_gdf_list[row_number]}')
    term_gdf.set_index('index',inplace=True)
    term_gdf.index.names = ['date']
    
    day = pd.Timedelta(1,'D')
    mid_date = date2 - ((date2 - date1)/2)
    day_multiplier = 7
    min_date = mid_date - (day_multiplier*day)
    max_date = mid_date + (day_multiplier*day)       
    
    x_atc = vel_pkl_full['x_atc']
    vel_pkl = vel_pkl_full.iloc[:,:-2]
    vel_sample_pkl = vel_pkl.loc[:,min_date:max_date] # velocity
    vel_dates = vel_sample_pkl.loc[:,~vel_sample_pkl.T.duplicated(keep='first')]
    
    vel_term_sample = []
    for col in vel_dates:
        if col in term_gdf.index:
            term_pt_vel = term_gdf.loc[col,:]
            x_pt_vel = list(term_pt_vel.geometry.coords)[0][0]
            # print(x_pt_vel)
            
            term_idx_vel = x_atc[x_atc == x_pt_vel].index[0]
            vel_sample = vel_dates.loc[term_idx_vel+300:term_idx_vel+900,col].mean()
            vel_term_sample.append(vel_sample)
        else: 
          vel_term_sample.append(np.nan)    
    
    days_from_calving = pd.Series(vel_dates.columns) - mid_date
    
    
    date_from_ax2 = mdates.DateFormatter('%m-%d')
    # print(vel_term_sample)
    percent_change = (vel_term_sample - np.nanmean(vel_term_sample))*100
    vel_sample_data_days = days_from_calving + pd.Timedelta(12,'h') # added 12 to get mid date for 24 hour vels #/day 
    ax.plot(vel_sample_data_days/day,vel_term_sample,'-', markeredgecolor='k',c='tab:gray',
            alpha=0.2,zorder=1,label='Terminus velocity')
    # ax.set_xlim(min_date,max_date)
    # ax.set_xticklabels(ax[2].get_xticks(), rotation = 30, ha="right")
    # ax.xaxis.set_major_formatter(date_from_ax2)
    ax.set_ylim(15,38)
    ax.grid(linestyle='--')
        # ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    ax.set_ylabel('m d$^{-1}$',fontsize=fontsize)
    # ax.set_ylabel('percent change',fontsize=fontsize)
    ax.set_xlabel('Days around calving event',fontsize=fontsize)
    
    # average_vel_arr[:len(vel_term_sample),i] = vel_term_sample
    # average_day_arr[:len(vel_sample_data_days),i] = vel_sample_data_days
    i+=1
    days_from_calving_df = pd.DataFrame(vel_sample_data_days)
    days_from_calving_df[date1] = vel_term_sample
    # date_1_fn = str(date1).replace(' ','_').replace(':','_')
    # date_2_fn = str(date2).replace(' ','_').replace(':','_')
    
    df_store.append(days_from_calving_df)
    
# calving_signiture = np.nanmean(average_vel_arr,axis=0)
# calving_signiture_days = np.nanmean(average_day_arr,axis=0)

df_all = pd.concat(df_store) # combine all the velocity sample dataframes
df_all.rename(columns={0:'days_around_calving'},inplace=True)
df_all['days'] = df_all['days_around_calving'].dt.days

pkl_vel_out = '/media/laserglaciers/upernavik/flexure_manuscipt/pkls/'



df_all.rename(columns={0:'days_around_calving'},inplace=True)
df_all['days'] = df_all['days_around_calving'].dt.days
df_all.to_pickle(f'{pkl_vel_out}vel_non_flexure_style_calving.pkl')

# gb = df_all.groupby(['days'],dropna=True).mean()
# calving_signiture = np.nanmean(gb.iloc[:,1:],axis=1) # I removed the days_around_calving column

# ax.plot(gb.index.values,calving_signiture,'-', markeredgecolor='k', c='k',
#         alpha=1,zorder=2,label='Terminus velocity')


# ax.tick_params(axis='both',which='major',labelsize=20)
# ax.set_title(f'Non-Flexure Calving Signiture', size=20)
# file_name = f'nonflexure_calving_signiture'
    
# plt.savefig(f'{op}{file_name}.png',dpi=300)
    
    
    
    # plt.savefig(f'{op}{file_name}.png')
    
        # ax.set_title('Velocity Sample')
        # ax.xaxis.set_minor_locator(MultipleLocator(1))
        
    
    
    
    
    
    
    
    
    
    