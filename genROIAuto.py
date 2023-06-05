#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   genROIAuto.py
@Time    :   2022/11/02 18:45:50
@Author  :   nico 
@Version :   1.0
@Description: 

This code will automatically take generated ribbonsn from Matlab in their proper directory and generate 
ROIs that can be used with the PWS analysis software. Instructions to run:

1. conda activate into your pwpsy environment 
2. run python code using 
    - python genROIAuto.py 
3. Follow prompt to select the folder where your Cell data is.

for any questions please reach out to Nico!
'''

import os 
import matplotlib.pyplot as mpl
import matplotlib.pyplot as plt
import cv2 as cv
import numpy as np
import pandas as pd
import h5py
from tkinter import Tk
from tkinter.filedialog import askdirectory 
import json
import tifffile as tif 
import pwspy.dataTypes as pwsdt
import rasterio 


### FUNCTIONS ###
def importMask(maskArray,cellFolder,roiName,roiNumber):
    """
    This function will generete a mask to be used with a particular 
    acquisiton object. 

    Args:
        maskArray (nd array ): array of mask, must be of type rasterio.bool_
        cellFolder (str): folder where the acquisition is
        roiName (int): integer for roi name ( if number if not can use string name as well)
    """
    roi = pwsdt.Roi.fromMask(maskArray)
    acq = pwsdt.Acquisition(cellFolder)
    acq.saveRoi(f'{roiName}',roiNumber,roi)
    print(f"uploaded {roiName} {roiNumber}")
    return 

def generateROIs(cellPath):
    ''' 
    This funciton will take in the path for a specific cell folder
    and generate rois for all of the tif files for the ribbons generated
    by the Matlab code 
    '''
    tifList =  [ file for file in os.listdir(cellPath) if file.endswith('.tif')] 
    #print(tifList)
    acquisitionPath = cellPath
    for i in range(len(tifList)):
        # this part generates the roi name that will be used for the descriptor and the roi number to save that roi ribbon as a number
        t = tifList[i].split('_')
        #print(t)
        #getting ROI Number right 
        #newROINum =  t[3][:-1]+str(int(t[3][-1])-1)
       # print(newROINum)
        #print(i)
        if t[1] == 'nuc' or  t[1] == 'nof':
            roiName = f'{t[1]}_ribbons_{t[3]}'
            roiNumber = i - 7*int(i/7)+1
            #print(roiName,roiNumber)
            # now we will read in the images at that path and and rotate 
            tifPath = os.path.join(cellPath,tifList[i])
            arrayImage = tif.imread(tifPath)
            maskArray = np.copy(arrayImage).astype(rasterio.bool_).T # T left in because dont have updated code that rotates the images in matlab, will remove later 
            # call import mask to create mask for all the masks 
            try:
                importMask(maskArray,acquisitionPath,roiName,roiNumber)
                pass
            except IndexError:
                print(f'issue with index in from mask for {roiName}')
                i +=1 
                continue

    tifList = []
    #print(tifList)
    return 

dataPath = askdirectory(title='Select Data folder where PWS Cell data is') # shows dialog box and return the path

## have to sort the folders in order then remove any other .DS_store files 
fileList = [] 
for file in os.listdir(dataPath):
    if  not len(file) == 7 and file[4] != '9':
        if file == ".DS_Store":
            # removing DS.STORE 
            os.remove(os.path.join(dataPath,file))
        else: 
            fileList.append(file)
        cellPath = os.path.join(dataPath,file)
        print(cellPath)
        #print(os.listdir(cellPath))
        try:
            generateROIs(cellPath)
            pass
        except OSError:
            print(f'file already exists {file} ')
            continue
    