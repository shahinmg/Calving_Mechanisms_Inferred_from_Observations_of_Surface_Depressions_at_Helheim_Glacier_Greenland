#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 18:24:29 2024

@author: laserglaciers
"""

import numpy as np
import os
import rasterio
import scipy
import scipy.ndimage as ndimage
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel




# https://blog.finxter.com/effective-gaussian-filtering-of-images-with-nans-in-python-using-matplotlib/
def nan_gaussian_filter(image, sigma=1):
    V = np.where(np.isnan(image), 0, image)  # Values (with NaNs to 0)
    V[np.isnan(image)] = 0  # Ignore NaNs
    VV = ndimage.gaussian_filter(V, sigma=sigma)
    W = np.where(np.isnan(image), 0, 1)  # Weights
    WW = ndimage.gaussian_filter(W, sigma=sigma)
    return VV/WW



def strain_rates(Ux, Vy, spacing, smooth=False, smooth_method = 'scipy', x_std = 1, y_std = 1,
                 smooth_kernal=1, rotate=False):
    
    dx, dy = spacing
    
    Ux_grad = np.gradient(Ux)
    Vy_grad = np.gradient(Vy)
    Exx = Ux_grad[1]/(dx)
    Eyy = (Vy_grad[0]/(dy))
    Exy = (0.5*((Ux_grad[0]/dy)+(Vy_grad[1]/dx)))
    
    
    Exx = np.ma.masked_array(Exx, np.isnan(Exx))
    Eyy = np.ma.masked_array(Eyy, np.isnan(Eyy))
    Exy = np.ma.masked_array(Exy, np.isnan(Exy))

    
    if smooth:
        
        if smooth_method == 'scipy':
            
            Exx = ndimage.gaussian_filter(Exx, sigma=smooth_kernal)
            Eyy = ndimage.gaussian_filter(Eyy, sigma=smooth_kernal)
            Exy = ndimage.gaussian_filter(Exy, sigma=smooth_kernal)
            
        elif smooth_method == 'astropy':
            kernel = Gaussian2DKernel(x_stddev = x_std, y_stddev = y_std, 
                                      x_size = smooth_kernal, y_size = smooth_kernal)
            
            
            Exx = convolve(Exx, kernel) #smooth it with Gaussian
            Eyy = convolve(Eyy, kernel)
            Exy = convolve(Exy, kernel)
            
    
    
    if rotate:
        
        theta = np.arctan2(Vy,Ux) 

        #---- Rotate Strain Tensors from pg 352 vdv 2nd edition          
        Exx = (Exx*(np.cos(theta)**2)) + (Eyy*((np.sin(theta)**2))) + (2*Exy*(np.sin(theta))*(np.cos(theta)))
        Eyy = (Exx*(np.sin(theta)**2)) + (Eyy*((np.cos(theta)**2))) - (2*Exy*(np.sin(theta))*(np.cos(theta)))
        Exy = (-Exx*np.sin(theta)*np.cos(theta)) + (Eyy*np.sin(theta)*np.cos(theta))+(Exy*(2*(np.cos(theta)**2) - (1))) # apply double-angle formula to cos**2 - sin**2
        # why is exx negative here? 
        
        strain_dict = {'Exx': Exx,
                       'Eyy': Eyy,
                       'Exy': Exy,
                       'theta': theta
                        }
        
        return strain_dict
        
    strain_dict = {'Exx': Exx,
                   'Eyy': Eyy,
                   'Exy': Exy
                   }
                        
    return strain_dict


def resistive_stress(Exx, Eyy, Exy, surf, bed, spacing, theta, n=3, B=500, BL=500,
                     smooth=False, smooth_surf=False, smooth_method = 'scipy', 
                     smooth_kernal=1, scipy_mode = 'reflect', x_std = 1, y_std = 1):
    """
    Compute strain fields and return in dictionaries

    Parameters
    ----------
    Exx : ndarray
        2D array of extensional strain rates. Most likely from the strain_rates function. 
    Eyy : ndarray
        2D array of yy strain rates. Most likely from the strain_rates function.
    Exy : ndarray
        2D array of shear strain rates. Most likely from the strain_rates function.
    surf : ndarray
        DEM surface.
    bed : ndarary
        Bed arrary.
    spacing : tuple
        tuple of the spacing used for the gradient. Try to get from the velocity raster's transformation
    theta : float
        Flow angle direction. Most likely from the strain_rates function.
    n : int, optional
        Glen's flow law exponent. The default is 3.
    B : float, optional
        Glen's flow law rate parameter in units of kPa yr^1/3. The default is 400.
    BL : TYPE, optional
        DESCRIPTION. The default is 300.
    smooth : bool, optional
        Smooth data using scipy.ndimage.gaussian_filter. The default is False.
    smooth_surf : bool, optional
        Smooth surface data using scipy.ndimage.gaussian_filter. The default is False.
    smooth_kernal : int, optional
        kernal size to use for scipy.ndimage.gaussian_filter. The default is 1.

    Returns
    -------
    resistive_stress_dict : dict
        DESCRIPTION.

    """
    
    ## ----------------------------
    ## CALCULATE RESISTIVE STRESSES
    ## ----------------------------

    # Need to rename these variables 
    
    if smooth_surf:
        surf = ndimage.gaussian_filter(surf, 2, mode=scipy_mode)
    
    H = surf - bed # get ice thickness
    H = np.ma.masked_array(H, np.isnan(H))
    theta = np.ma.masked_array(theta, np.isnan(theta))
    
    Ezz = -Exx - Eyy
    Eeff = pow(.5*(pow(Exx,2)+pow(Eyy,2)+pow(Ezz,2))+pow(Exy,2),.5)
    Eeff23 = pow(Eeff,-.666)
    
    Oxx = B * Exx * Eeff23
    Oyy = B * Eyy * Eeff23
    Oxy = BL * Exy * Eeff23 # change Bl to equal B
    
    Oxy_rot = (-Oxx*np.sin(theta)*np.cos(theta)) + (Oyy*np.sin(theta)*np.cos(theta))+(Oxy*(2*(np.cos(theta)**2) - (1)))
    
    Rxx = (2 * Oxx) + Oyy
    Ryy = (2 * Oyy) + Oxx
    Rxy = Oxy
    
    HRxx = Rxx * H
    HRxy = Rxy * H
    HRyy = Ryy * H
    HRxy_rot = (-HRxx*np.sin(theta)*np.cos(theta)) + (HRyy*np.sin(theta)*np.cos(theta))+(HRxy*(2*(np.cos(theta)**2) - (1)))

    dx,dy = spacing
    ## ----------------------------
    ## CALCULATE STRESS GRADIENTS 
    ## ----------------------------
    HRxxgrad = np.gradient(HRxx)
    HRxygrad = np.gradient(HRxy)
    HRyygrad = np.gradient(HRyy)
    
    ddxHRxx = HRxxgrad[1]/(dx)
    ddyHRxx = HRxxgrad[0]/(dy)
    ddxHRxy = HRxygrad[1]/(dx)
    ddyHRxy = HRxygrad[0]/(dy)
    ddxHRyy = HRyygrad[1]/(dx)
    ddyHRyy = HRyygrad[0]/(dy)

    # Resistance is +  
    # F_long_along = (np.cos(theta)**3)*(ddxHRxx) + ((np.cos(theta)*(np.sin(theta)**2))*ddxHRyy) + (2*(np.cos(theta)**2)*np.sin(theta)*ddxHRxy) + \
    # (((np.cos(theta)**2)*np.sin(theta))*ddyHRxx) + ((np.sin(theta)**3)*ddyHRyy) + (2*np.cos(theta)*(np.sin(theta)**2)*ddyHRxy)
    
    # F_long_across = (-np.sin(theta)**3)*(ddxHRxx) - (((np.cos(theta)**2)*np.sin(theta))*ddxHRyy) + (2*np.cos(theta)*(np.sin(theta)**2)*ddxHRxy) + \
    # (np.cos(theta)*((np.sin(theta)**2))*ddyHRxx) + ((np.cos(theta)**3)*ddyHRyy) - (2*(np.cos(theta)**2)*np.sin(theta)*ddyHRxy)
    
    # F_lat_along = (-(np.cos(theta)*(np.sin(theta)**2))*(ddxHRyy - ddxHRxx)) - (np.sin(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddxHRxy) + \
    # ((np.cos(theta)**2)*(np.sin(theta))*(ddyHRyy-ddyHRxx)) + (np.cos(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddyHRxy)
    
    # F_lat_across = ((np.cos(theta)**2)*np.sin(theta)*(ddxHRyy - ddxHRxx)) + (np.cos(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddxHRxy) + \
    # (np.cos(theta)*(np.sin(theta)**2)*(ddyHRyy-ddyHRxx)) + (np.sin(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddyHRxy)
    
    F_long_along = (np.cos(theta)**3)*(ddxHRxx) + (((np.sin(theta)**2)*np.cos(theta))*ddxHRyy) \
    + (2*np.sin(theta)*(np.cos(theta)**2)*ddxHRxy) + (((np.cos(theta)**2)*np.sin(theta))*ddyHRxx) \
    + ((np.sin(theta)**3)*ddyHRyy) + (2*np.cos(theta)*(np.sin(theta)**2)*ddyHRxy)
    
    F_long_across = (-np.sin(theta)**3)*(ddxHRxx) - (((np.cos(theta)**2)*np.sin(theta))*ddxHRyy) \
    + (2*np.cos(theta)*(np.sin(theta)**2)*ddxHRxy) + (np.cos(theta)*((np.sin(theta)**2))*ddyHRxx) \
    + ((np.cos(theta)**3)*ddyHRyy) - (2*(np.cos(theta)**2)*np.sin(theta)*ddyHRxy)
    
    F_lat_along = (-(np.cos(theta)*(np.sin(theta)**2))*(ddxHRyy - ddxHRxx)) - \
    (np.sin(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddxHRxy) + \
    ((np.cos(theta)**2)*(np.sin(theta))*(ddyHRyy-ddyHRxx)) + \
    (np.cos(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddyHRxy)
    
    F_lat_across = ((np.cos(theta)**2)*np.sin(theta)*(ddxHRyy - ddxHRxx)) + \
    (np.cos(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddxHRxy) + \
    (np.cos(theta)*(np.sin(theta)**2)*(ddyHRyy-ddyHRxx)) + \
    (np.sin(theta)*((np.cos(theta)**2)-(np.sin(theta)**2))*ddyHRxy)
    

    
    if smooth:
        
        if smooth_method == 'scipy':
            
            F_long_along = ndimage.gaussian_filter(F_long_along, sigma=smooth_kernal)
            F_long_across = ndimage.gaussian_filter(F_long_across, sigma=smooth_kernal)
            F_lat_across = ndimage.gaussian_filter(F_lat_across, sigma=smooth_kernal)
            F_lat_along = ndimage.gaussian_filter(F_lat_along, sigma=smooth_kernal)

            
        elif smooth_method == 'astropy':
            kernel = Gaussian2DKernel(x_stddev = x_std, y_stddev = y_std, 
                                      x_size = smooth_kernal, y_size = smooth_kernal)
            
            
            F_long_along = convolve(F_long_along, kernel)
            F_long_across = convolve(F_long_across, kernel)
            F_lat_across = convolve(F_lat_across, kernel)
            F_lat_along = convolve(F_lat_along, kernel)
            
            

    slope = np.gradient(surf)
    slope_x = slope[1]/dx
    slope_y = slope[0]/dy
    slope_along = (slope_x*np.cos(theta)) + (slope_y*(np.sin(theta)))
    slope_across = (slope_y*np.cos(theta)) - (slope_x*(np.sin(theta)))
    
    
    
    # Need to conver to negative to make resistance postive  
    resistive_stress_dict = {'F_long_along': -F_long_along,
                             'F_long_across': -F_long_across,
                             'F_lat_across': -F_lat_across,
                             'F_lat_along': -F_lat_along,
                             'H': H,
                             'slope_x': slope_x,
                             'slope_y': slope_y,
                             'Rxx': Rxx,
                             'Ryy': Ryy,
                             'Rxy': Rxy,
                             'surface':surf
                             }
    
    return resistive_stress_dict



def driving_stress(F_long_along, F_lat_along, F_lat_across, F_long_across, H, slope_x, slope_y,
                   theta, rho=0.917, rhow=1.0, g=9.8, 
                   smooth=False, smooth_method = 'scipy', smooth_kernal=1, scipy_mode = 'reflect',
                   x_std = 1, y_std = 1):
    
    
    H = np.ma.masked_array(H, np.isnan(H))
    theta = np.ma.masked_array(theta, np.isnan(theta))
    
    slope_x = np.ma.masked_array(slope_x, np.isnan(slope_x))
    slope_y = np.ma.masked_array(slope_y, np.isnan(slope_y))
    
    Tdx = -rho * g * H * slope_x
    Tdy = -rho * g * H * slope_y
    
    Td_along = Tdx*(np.cos(theta)) + Tdy*(np.sin(theta))     
    Td_across = Tdy*(np.cos(theta)) - Tdx*(np.sin(theta))   
    
    if smooth:
        
        if smooth_method == 'scipy':
            
            Td_along = ndimage.gaussian_filter(Td_along, sigma=smooth_kernal) # I smooth it here to avoid the ribbing effects
            Td_across = ndimage.gaussian_filter(Td_across, sigma=smooth_kernal) # that are described by Olga (2013). 9?


        elif smooth_method == 'astropy':
            
            kernel = Gaussian2DKernel(x_stddev = x_std, y_stddev = y_std, 
                                      x_size = smooth_kernal, y_size = smooth_kernal)
            Td_along = convolve(Td_along, kernel)
            Td_across = convolve(Td_across, kernel)  
    
    Td_along =  np.ma.masked_array(Td_along, np.isnan(Td_along))
    F_long_along =  np.ma.masked_array(F_long_along, np.isnan(F_long_along))
    F_lat_across =  np.ma.masked_array(F_lat_across, np.isnan(F_lat_across))

    
    Tb_along = Td_along - F_long_along - F_lat_along # make sure all inputs are smoothed
    Tb_across = Td_across - F_long_across - F_lat_across
    
    # if smooth:
    #     if smooth_method == 'scipy':
            
    #         Tb_along = nan_gaussian_filter(Tb_along, sigma=smooth_kernal) 
    #         Tb_across = nan_gaussian_filter(Tb_across, sigma=smooth_kernal)

    
    #     elif smooth_method == 'astropy':
    #         kernel = Gaussian2DKernel(x_stddev = 1, y_stddev = 1, 
    #                                   x_size = smooth_kernal, y_size = smooth_kernal)
    #         Tb_along = convolve(Tb_along, kernel)
    #         Tb_across = convolve(Tb_across, kernel)  
    
    # rho = .917 #kg/m2
    # rhoS = 1.02
    # Hab = thick + (rhoS/rho)* bed #height above buoyancy
    # Udef = 0.5* thick*(Tb_along / B)**3 #deformational velocity
    
    driving_stress_dict = {'Tdx': Tdx,
                           'Tdy': Tdy,
                           'Td_along': Td_along,
                           'Td_across': Td_across,
                           'Tb_along': Tb_along,
                           'Tb_across': Tb_across
                            }
    
    return driving_stress_dict