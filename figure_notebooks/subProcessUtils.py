#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:26:45 2020

@author: laserglaciers
"""
import subprocess
import os
import numpy as np



def date_sort(path,extension,add_20=True):
    """
    Parse through atlas laz filenames to create np.datetime64 from filename 
    objects and make list of tuples with original file name and datetime.
    
    Parameters
    ----------
    path: String
        Path to where laz/tif files are saved
    
    extension: String
        File extension for specific file type to be listed
    
    Returns
    -------
    dates_sort\n
    A sorted list of tuples with laz/tif file in first position and np.datetime64
    in the second posistion. List sorted in ascending order.\n
    fileList_sort
    sorted file list\n
    npDates_sort
    sorted np dates
    """
    fileList = [f for f in os.listdir(path) if f.endswith(extension)]
    
    if add_20 == True:
        dates = ['20'+date.replace('_','').replace('-','')[:-len(extension)] for date in fileList] #remove '_" and add 20
        
    else:
        dates = [date.replace('_','').replace('-','')[:-len(extension)] for date in fileList] #remove '_" and add 20
        
    npDates = [np.datetime64(date[:4]+'-'+date[4:6]+'-'+date[6:8]+' '+
            date[8:10]+':'+date[10:12]+':'+date[12:14]) for date in dates]
    
    
    dates_sort = sorted(list(zip(fileList,npDates)),key = lambda t:t[1])
    fileList_sort, npDates_sort = zip(*(dates_sort))
    
    return dates_sort, fileList_sort, npDates_sort



def date_sort_flexures(path,extension):
    """
    Parse through atlas laz filenames to create np.datetime64 from filename 
    objects and make list of tuples with original file name and datetime.
    
    Parameters
    ----------
    path: String
        Path to where laz/tif files are saved
    
    extension: String
        File extension for specific file type to be listed
    
    Returns
    -------
    dates_sort\n
    A sorted list of tuples with laz/tif file in first position and np.datetime64
    in the second posistion. List sorted in ascending order.\n
    fileList_sort
    sorted file list\n
    npDates_sort
    sorted np dates
    """
    fileList = [f for f in os.listdir(path) if f.endswith(extension)]
    dates = [date.replace('-','').replace('_','')[:-len(extension)] for date in fileList] #remove '_" and add 20
    
    npDates = [np.datetime64(date[:4]+'-'+date[4:6]+'-'+date[6:8]+' '+
            date[8:10]+':'+date[10:12]+':'+date[12:14]) for date in dates]
    
    
    dates_sort = sorted(list(zip(fileList,npDates)),key = lambda t:t[1])
    fileList_sort, npDates_sort = zip(*(dates_sort))
    
    return dates_sort, fileList_sort, npDates_sort


def time_match(date_sort_zip,match_times):
    
    matched_files = []
    for tup in date_sort_zip:
        file,time = tup 
        
        if time not in match_times:
            pass
        else:
            matched_files.append(tup)
            
    return matched_files


def batch_process(laz_sort,svg_path,atlas_path,laz_path):
    """
    Runs ./atlas command for a directory of laz files and saves a displacement
    laz and svg file of vectors.
    
    Parameters
    ----------
    laz_sort : List
        A list of of zipped laz files at [0] and corresponding np.datetime64 
        object at [1].
    svg_path : String
        Path to where to store svg file. Laz path is specified in ATLAS.cpp
        
    atlas_path : String
        Path to where ./atals is built
    Returns
    -------
    None.

    """
    
    os.chdir(atlas_path)
    
    for i,laz in enumerate(laz_sort):
    
        if i == len(laz_sort)-1:
            break
        else:
            print('Calculating displacements between '+' '.join([laz[0],laz_sort[i+1][0]]))
            args = ['./atlas', laz_path+laz[0], laz_path+laz_sort[i+1][0], 
                    '--xshift=6.4 ', '--yshift=-.51' ]
            runAtlas = subprocess.Popen(args,stdout=subprocess.PIPE)
            runAtlas.wait()
            move_svg = subprocess.check_output(['mv', 'out.svg', svg_path+laz[0][:-4]+'.svg'], stdin=runAtlas.stdout)
        
    return None

