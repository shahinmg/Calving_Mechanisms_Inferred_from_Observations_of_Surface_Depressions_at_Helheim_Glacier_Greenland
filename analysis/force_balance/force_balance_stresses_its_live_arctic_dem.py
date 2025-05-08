#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:49:15 2024

@author: laserglaciers
"""

import os
import numpy as np
import scipy 
import force_balance as forces
import xarray as xr
import rasterio
import pandas as pd
from subProcessUtils import date_sort
from multiprocessing import pool, cpu_count


def nan_2_zero(array):
    array = np.where(np.isnan(array),0,array)
    return array

def force_balance_process(dem, bed, vel, dem_path, vel_path, op, smooth_kernal_size=3, 
                          B=500, BL=500, smooth_strain = False,  smooth_ressitive=False, smooth_ressitive_kernal_size=3,
                          smooth_driving=False, smooth_driving_kernal=3, smooth_method = 'scipy', x_std=1, y_std=1):
    
    
    dem_path = f'../resampled_vels_dems/arctic_dem_{resolution}m_crop/review/{dem}'
    vel_path = f'../resampled_vels_dems/selected_vels_{resolution}m/review/all_bands_zero_no_data_crop/{vel}'

    
    dem = xr.open_dataset(dem_path, engine='rasterio')
    bed = xr.open_dataset(bed_path, engine='rasterio')
    vel = xr.open_dataset(vel_path, engine='rasterio')
    # vel = vel.rio.interpolate_na(method = 'cubic')
    
    vx = vel.band_data.to_numpy()[1,:,:]
    vy = vel.band_data.to_numpy()[2,:,:]
    bed = bed.bed.to_numpy()[0,:,:]
    dem = dem.band_data.to_numpy()[0,:,:]
    
    # try turning nans to zero
    vx = nan_2_zero(vx)
    vy = nan_2_zero(vy)
    bed = nan_2_zero(bed)
    dem = nan_2_zero(dem)
    
    spacing = vel.rio.resolution()
    
    strain_dict = forces.strain_rates(vx, vy, spacing, smooth=smooth_strain, 
                                      smooth_kernal=smooth_kernal_size, rotate=True, smooth_method = smooth_method,
                                      x_std = x_std, y_std=y_std)
    Exx = strain_dict['Exx']
    Eyy = strain_dict['Eyy']
    Exy = strain_dict['Exy']
    theta = strain_dict['theta']
    
    
    
    resistive_stress_dict = forces.resistive_stress(Exx, Eyy, Exy, dem, bed, spacing, theta,
                                                    smooth=smooth_ressitive, smooth_kernal=smooth_ressitive_kernal_size, 
                                                    B=B, BL=BL,  smooth_method = smooth_method, x_std = x_std, y_std=y_std,
                                                    smooth_surf=True)
    
    
    F_long_along = resistive_stress_dict['F_long_along']
    F_lat_along = resistive_stress_dict['F_lat_along']
    F_lat_across = resistive_stress_dict['F_lat_across']
    F_long_across = resistive_stress_dict['F_long_across']
    H = resistive_stress_dict['H']
    slope_x = resistive_stress_dict['slope_x']
    slope_y = resistive_stress_dict['slope_y']
    Rxx = resistive_stress_dict['Rxx']
    Ryy = resistive_stress_dict['Ryy']
    Rxy = resistive_stress_dict['Rxy']
    surf = resistive_stress_dict['surface']
    
    driving_stress_dict = forces.driving_stress(F_long_along, F_lat_along, F_lat_across, F_long_across, 
                                                H, slope_x, slope_y, theta,smooth=smooth_driving, 
                                                smooth_kernal=smooth_driving_kernal,  smooth_method = smooth_method,
                                                x_std = x_std, y_std=y_std)
    
    Tdx = driving_stress_dict['Tdx']
    Tdy = driving_stress_dict['Tdy']
    Td_along = driving_stress_dict['Td_along']
    Td_across = driving_stress_dict['Td_across']
    Tb_along = driving_stress_dict['Tb_along']
    Tb_across = driving_stress_dict['Tb_across']
    
    balance = Td_along - F_long_along - F_lat_along - Tb_along
    # balance = tb - td + flong + flat
    
    force_balance_dict = {'Exx': Exx,
                          'Eyy': Eyy,
                          'Exy': Exy,
                          'theta': theta,
                          'F_long_along': F_long_along,
                          'F_lat_along': F_lat_along,
                          'F_lat_across': F_lat_across,
                          'F_long_across': F_long_across,
                          'Tdx': Tdx,
                          'Tdy': Tdy,
                          'Td_along': Td_along,
                          'Td_across': Td_across,
                          'Tb_along': Tb_along,
                          'Tb_across': Tb_across,
                          'slope_x': slope_x,
                          'slope_y': slope_y,
                          'Rxx': Rxx,
                          'Ryy': Ryy,
                          'Rxy': Rxy,
                          'surface': surf,
                          'H': H,
                          'balance':balance
                          }
    
    
    meta = {'width':vel.rio.width,
            'height': vel.rio.height,
            'count': len(force_balance_dict),
            'nodata': 0, #vel.rio.nodata need to encode this. I inspected vx manually and found nans
            'transform': vel.rio.transform(),
            'crs':vel.rio.crs,
            }
    
    
    with rasterio.open(f'{op}',mode='w',dtype=np.float32,**meta) as dst:
        
        for band_num, band_name in enumerate(force_balance_dict.keys(), start=1):
            band_data = force_balance_dict[band_name]
            dst.write(band_data,band_num)
            dst.set_band_description(band_num, band_name)
            
            

        
    
    print(f'saved {op}')
        

def date_match(vel_list, dem_list):
    matched_files = []
    
    for vel in vel_list:
        vel_year = pd.to_datetime(f'20{vel[:13]}'.replace('_',' '))
        
        for dem in dem_list:
            dem_year = pd.to_datetime(f'20{dem[:13]}'.replace('_',' '))
            
            if vel_year == dem_year:
                matched_files.append((vel,dem))
    
    return matched_files





smooth_kernal_size = 2
smoth_ressitive_kernal_size = 1.6
smooth_driving_kernal = 1.6
x_std = 7
y_std = 7

smooth_strain = True
smoth_ressitive = True
smooth_driving = True

# dem_smooth = 3
# vel_smooth = 8
B = 250
Bl = 210
resolution = 250

bed_path = f'../resampled_vels_dems/bedmachineV5_crop/review/bed_resampled_{resolution}m_zero_no_data.tif'
dem_path = f'../resampled_vels_dems/arctic_dem_{resolution}m_crop/review/'
vel_path = f'../resampled_vels_dems/selected_vels_{resolution}m/review/all_bands_zero_no_data_crop/'


op_root = f'/media/laserglaciers/upernavik/force_balances/its_live_arctic_dem_forces/forces_outfile/review/final/update_testing/'
op = f'{op_root}forces_out_smooth_strain_{smooth_strain}_{smooth_kernal_size}_smooth_ressitive_{smoth_ressitive}_{smoth_ressitive_kernal_size}_smooth_driving_{smooth_driving}_{smooth_driving_kernal}_nan_mask_smooth_surf_H_output_surface_rate_factor_{B}_Bl_{Bl}_CROP_proper_mask_no_vel_interp_v2/'


if not os.path.exists(op):
   os.makedirs(op)

def date_math(vel_list, dem_list):
    matched_files = []
    
    for vel in vel_list:
        vel_year = pd.to_datetime(vel[11:26]).year
        
        for dem in dem_list:
            dem_year = pd.to_datetime(dem[18:26]).year
            
            if vel_year == dem_year:
                matched_files.append((vel,dem))
    
    return matched_files

def its_live_s2_dates(file_name, kernal_size=0, no_smooth=True):
    date1 = file_name[11:26]
    date2 = file_name[108:123]
    
    out_file = f'{date1}_{date2}_force_balance_smooth_{kernal_size}.tif'
    
    return out_file
    


# vel_list = [vel for vel in os.listdir(vel_path) if vel.endswith('nc')]
vel_list = [vel for vel in os.listdir(vel_path) if vel.endswith('tif')]
dem_list = [dem for dem in os.listdir(dem_path) if dem.endswith('tif')]
matched_files = date_math(vel_list, dem_list)

for vel, dem in matched_files:
    
    out_file = its_live_s2_dates(vel,kernal_size=smooth_kernal_size)
    out_path = f'{op}{out_file}'
    force_balance_process(dem, bed_path, vel, dem_path, vel_path, out_path, smooth_kernal_size=smooth_kernal_size, 
                          B=B, BL=Bl, smooth_strain = smooth_strain,  smooth_ressitive=smoth_ressitive,
                          smooth_ressitive_kernal_size=smoth_ressitive_kernal_size, 
                          smooth_driving=smooth_driving, smooth_driving_kernal=smooth_driving_kernal, 
                          smooth_method = 'scipy', x_std = x_std, y_std = y_std)

