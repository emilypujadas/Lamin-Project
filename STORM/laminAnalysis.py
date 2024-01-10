#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   laminAnalysis.py
@Time    :   2024/01/09 16:22:17
@Author  :   Backman Lab - Nicolas Acosta 
@Version :   1.0
@Descr   :   Code should analyze lamin data and compare regional periphery 
             regional interior intensity signals within collected dataset.
             Code is tailored to internal lab organization of files, but can 
             be modified to work with any file organization schema.  
'''

# import functions 
# importing all necessary functions 
# importing functions 

# import functions 
# importing all necessary functions 
# importing functions 

import os 
import scipy.io as sio 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
from tkinter import Tk
from tkinter.filedialog import askdirectory 
import pandas as pd
#import h5py
import math
from mpl_toolkits import mplot3d
from matplotlib.animation import FuncAnimation

matplotlib.rcParams["figure.dpi"] = 100

### import functions ###
#importing packages
import tifffile as tif
import cv2 as cv 
import seaborn as sns

from scipy.spatial import distance
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from scipy.signal import savgol_filter as sf
from scipy.signal import find_peaks
from scipy import stats
from sklearn.cluster import DBSCAN 
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances_argmin_min
from sklearn.metrics import pairwise_distances
import skimage.exposure


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=18) 


# set path for figures  and data 
figuresPath = r'''\\backmanlabnas.myqnapcloud.com\Public\Nico\Emily_Collab\LaminManuscript\figures'''
dataSavePath = r'''\\backmanlabnas.myqnapcloud.com\Public\Nico\Emily_Collab\LaminManuscript\processedData'''


# necessary functions 
def readImages(imgPath):
    ''' 
    read image file and return array, dtype and the  shape 
    '''
    img = cv.imread(imgPath,0) # read image as grayscale
    try: 
        imgShape = img.shape
    except AttributeError:
        print("Initial attempt with opencv failed, will try again with tifffile package...")

    if img is None:
        img2 = tif.imread(imgPath)
        img = img2 
        imgShape = img2.shape
        
    ## create return obj 
    outList= [img,type(img),imgShape]
    return outList


def generateMask(stormIm,pwsMap = None ,dfactor=18,efactor=15,kernel=np.ones((3,3),np.uint8),q=0.99,iter = 39,diff=5,sig=15):
    ''' 
    This function takes in the storm Image and generates the mask used
    to mask each image, need to play around with the dilation and erosion
    amounts to remove bacgkrounds. It uses dilation erosion and connected areas
    to remove small noise localizations. This is easier if labeling is good. 
    ''' 
    #### create the mask for storma nd mask the dmap and storm image that way #### 
    # have the storm and the dmap, no masks so code skips the dmap warp and the upsampled dmap steps 
    # mask it first with map large to get rid of the backgroudn 
    matplotlib.rcParams["figure.dpi"] = 100 # for plotting  
    #stormIm = stormIm[:,:stormIm.shape[0]-stormIm.shape[1]] # get it square or not 
    # take in loaded storm and mask it with larger fluorescence image if you have that or with the ROI drawn by PWS
    testStorm = cv.convertScaleAbs(stormIm,alpha=10)
    testStorm = testStorm.astype(np.uint8) 

    if type(pwsMap) != None:
        stormIm2 = cv.bitwise_and(testStorm,testStorm,mask=pwsMap).astype(np.uint8)
    else:
        stormIm2 = testStorm

    # create the mask and mask the dmap image and the storm image 
    # create mask from STORM reconstructed image 
    th2,threshSTORM = cv.threshold(stormIm2.astype(np.uint8),0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
    threshSTORM = threshSTORM.astype(np.uint8) 
    # blur the  image 
    stormMask = cv.GaussianBlur(threshSTORM,(5,5),0)

    #############################################################################################################################
    # dialte 
    # threshold the dmap image by generating STORMROI MASK BELOW
    kernel = np.ones((3,3),np.uint8)
    dstormMask = cv.dilate(src=stormMask,dst=None,kernel=kernel,iterations=dfactor)
    dStormMask = cv.erode(src=dstormMask,dst=None,kernel=kernel,iterations=efactor+5)# threshold the dmap image by generating STORMROI MASK BELOW

    ############################################################################################################################## do connectedCompononents  and filter 
    # perform conenctedComponents and filter the sizes to filter
    c = 8
    nb_blobs, im_with_separated_blobs, stats, _ = cv.connectedComponentsWithStats(dStormMask,c)
    #plt.imshow(dStormMask)
    #plt.show()
    #plt.imshow(im_with_separated_blobs)
    #plt.show()

    # filter through labels -- first check thresh
    sizes = stats[:,-1] 
    sizes = sizes[1:]
    # set placeholder to only keep blobs of interest
    stormMask = np.zeros_like(im_with_separated_blobs)
    min_size =  np.quantile(sizes,q) #np.mean(sizes) 
    # iterate through blobs and make sure blob size is large than min  
    for blob in range(len(sizes)):
        if sizes[blob]>=min_size:
            stormMask[im_with_separated_blobs == blob +1] = 255

    stormMask = stormMask.astype(np.uint8) 
    #plt.imshow(stormMask)
    #plt.show()

    # final dilation and eroisoon to fill the holes 
    iter = iter
    tDilate = cv.dilate(src=stormMask,dst=None,kernel=kernel,iterations=iter)
    tErode = cv.erode(src=tDilate,dst=None,kernel=kernel,iterations=iter-diff)
    #plt.imshow(tErode)
    #plt.show()

    blur = cv.GaussianBlur(tErode,(0,0),sigmaX = sig , borderType = cv.BORDER_DEFAULT )
    stormMask = skimage.exposure.rescale_intensity(blur, in_range=(120,255), out_range=(0,255))
    stormMask = stormMask.astype(np.uint8)

    #plt.imshow(stormMask),plt.axis('off')
    #plt.show()

    return stormMask,blur,tErode 



def padImages(allSlices,allImages,padDist=100):
    ''' 
    As func name states will pad the images by a certain amount such that 
    the dilation does not affect the masks --> its an attempt to avoid the limitations
    '''
    ## all slices code 
    allSlices2 = [] 
    allImages2 = []
    for i in range(len(allSlices)):
        # index data 
        tImg = allSlices[i] 
        # get shapes
        h,w = tImg.shape[:2]
        hNew,wNew = h+padDist*2,w+padDist*2
        padding_value = 0 

        ## create padded zeros array of desired size
        padData = np.zeros((hNew,wNew))

        # pad array by indexing the data as such 
        padData[padDist:h+padDist,padDist:w+padDist] = tImg

        allSlices2.append(padData)

        # all Images code 
        tImg2 = allImages[i]
        
        # transpose it to iterate in the way i like lol
        if tImg2.shape[0]<1000:
            tImg2 = tImg2.T
            #print(tImg2[:,:,0].shape)
        
        # we have to make a cube of zeros of desired shape. this just makes it easier when 
        # masking later as we can map the function to the cube list 
        tempCube = np.zeros((hNew,wNew,3))
        
        
        # iterate through the cube and pad, then save to the correspinding shape 
        for s in range(3):
            # get the slice in the image 
            tempSlice = tImg2[:,:,s]
            #print(tempSlice.shape)
            # use shape and new positions as the same from before 
            h,w = tImg.shape[:2]
            hNew,wNew = h+padDist*2,w+padDist*2
            # create new zeros array each time of size 
            padedData = np.zeros((hNew,wNew))
           # print(padedData.shape)
            # pad array by indexing the data as such 
            padedData[padDist:h+padDist,padDist:w+padDist] = tempSlice.T
            # now save tempslice to a slice in teh cube and append to the list of cubes 
            tempCube[:,:,s] = padedData 
        allImages2.append(tempCube.T)


    return allSlices2,allImages2

# make function for generating the masks and keeping them 
def genMasks(allSlices):
    '''
    take the fileList and then generate masks from the second channel image 
    '''
    # array for 
    maskArray = []
    # now do this for all iamges in the folder 
    for s in range(len(allSlices)):
        # now lets do a test thresholding and mask making so we can mask all the other channels in the image 
        t = allSlices[s]
        # threshold 
        th,threshImg =  cv.threshold(t.astype(np.uint8),0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
        # fill holes with dilation of mask 
        kernel = np.ones((3,3),np.uint8)
        # dilate/erosion
        dMask = cv.dilate(src=threshImg,dst=None,kernel=kernel,iterations=10)
        eMask = cv.erode(src=dMask,dst=None,kernel=kernel,iterations=10)
        # append to  maskArray
        maskArray.append(eMask)

        if s == 1 :
            # plot image 
            plotData = [threshImg,dMask,eMask]
            fig,axs = plt.subplots(nrows=1,ncols=3,figsize=(15,7))
            for i,ax in enumerate(axs.flatten()):
                ax.imshow(plotData[i]),ax.axis('off')
                if i == 0:
                    ax.set_title(f'Image {s} in directory')

    return maskArray


def genPeriphRegions(maskArray,iter =2):
    '''takes in maskArray and then outputs 
        a list of lists where each enetry is the 
        three masks of interest together 
    '''
    kernel = np.ones((3,3))
    regionalList = []
    for i in range(len(maskArray)):
        t = maskArray[i]    
        contT,_= cv.findContours(t,cv.RETR_EXTERNAL,cv.CHAIN_APPROX_NONE)
        contImg = np.zeros_like(t)
        #draw contours then display 
        contImg = cv.drawContours(contImg,contT,-1,(255),1)
        # dilate the  and then bitwise_and to only keep internal region of the image for calculations 
        kernel = np.ones((3,3),np.uint8) 
        dcontImg = cv.dilate(contImg,kernel,iterations=iter)
        # mask the dilated contour mask in order to keep whatever is inside and not what is outside 
        periphMask = cv.bitwise_and(dcontImg,t)
        innerRegion = t - periphMask
        # add to a a list ot keep track of it all 
        regionalList.append([t,periphMask,innerRegion])

    return regionalList

# now define mult function 
def multImg(mask,j,allSlices):
    ts = allSlices[j]
    return ts * mask


def calcIntensity(allSlices,regionalList,ct_idx = 5):
    ''' 
    Will go through each image after partitioned for treatment and report
    the proportion of k9 signal that lies within the periphery vs the
    interior region
    '''
    calcSum = lambda x: np.sum(x[x!=0])
    # iterate through each image, mask by each and then save that masked value for data keeping
    maskedAll = [ [multImg(regionalList[j][i],j,allSlices) for i in range(len(regionalList[j]))] for j in range(len(regionalList)) ]
    # make anonymous function to calculate the sum of nonzero values within the image 
    sumArray = [ list(map(calcSum,sublist)) for sublist in maskedAll ] 
    # now divide get resultArray by normalizing by total 
    resultsAll  = [list(map(lambda x: x/sublist[0],sublist)) for sublist in sumArray]

    ctrlRes = resultsAll[ct_idx:]
    treatRes = resultsAll[:ct_idx]

    return maskedAll,sumArray,[ctrlRes,treatRes]

def calcIntensity_v2(allSlices,regionalList,ct_idx = 5):
    ''' 
    Will go through each image after partitioned for treatment and report
    the proportion of k9 signal that lies within the periphery vs the
    interior region --> added normalization by count first before normalizing by the total nuclear intensity 

    so its sum(pixels)/size_mask / total nuclear sum / size mask 
    '''
    calcSum = lambda x: np.sum(x[x!=0])
    getSize = lambda x: len(x[x!=0])
    # iterate through each image, mask by each and then save that masked value for data keeping
    maskedAll = [ [multImg(regionalList[j][i],j,allSlices) for i in range(len(regionalList[j]))] for j in range(len(regionalList)) ]
    # make anonymous function to calculate the sum of nonzero values within the image 
    sumArray = [ list(map(calcSum,sublist)) for sublist in maskedAll ] 
    maskSizeArray  = [ list(map(getSize,sublist)) for sublist in maskedAll ] 

    # make life easier 
    normRes = []
    for i in range(len(maskedAll)): 
        normRes.append(list(map(lambda x,y:x/y,sumArray[i],maskSizeArray[i])))

    # now divide get resultArray by normalizing by total 
    resultsAll  = [list(map(lambda x: x/sublist[0],sublist)) for sublist in normRes]

    ctrlRes = resultsAll[ct_idx:]
    treatRes = resultsAll[:ct_idx]

    return maskedAll,sumArray,normRes,[ctrlRes,treatRes]


# load in the data 
dataPath = askdirectory()   
    
 
        
        
# get paths in a folder for processing 
fileList = [ file for file in os.listdir(dataPath) if file.endswith('.tif') ]
controlList = [file for file in fileList if file.startswith('Control')]
treatList =   [file for file in fileList if file.startswith('Auxin')]

fileList = [file for file in fileList if len(file)<24] 

# get only the first of each and pad it so we can do the mask genreation without any issues 
allImages = []
allSlices = [] 

for i in range(len(fileList)):
    path = os.path.join(dataPath,fileList[i])
    img,imgType,imgShape = readImages(path)
    allImages.append(img)
    allSlices.append(img[0])

# pad the images 
allSlices,allImages = padImages(allSlices,allImages)
maskArray = []
#gen massks 
for i in range(len(allSlices)):
    # gen mask and then binarize it 
    tempMask = generateMask(allSlices[i])[0]
    # binarize and add to list 
    t,tempBinary = cv.threshold(tempMask,0,255,cv.THRESH_BINARY)
    maskArray.append(tempBinary)
#plt.imshow(maskArray[0]),plt.show()

# Get Laminar Region
# now  make function to generate the laminarMask and store it
# within a list where you have nucMask,PeriphMask,innerMask
regionalList = genPeriphRegions(maskArray,iter=5) # starting with iter 10 to see how it affects results 

maskedAll,sumArray,normRes,resultsAll = calcIntensity_v2(allSlices,regionalList,ct_idx=5)

#  reformat data and plot results  # 
# plot as average barplot with error bars using seaborn 
# save data as csv dataframe such taht it can be used in graphpad or any other
# plotting software 

cR,tR = resultsAll   
ctrl_df = pd.DataFrame(cR,columns=['TotalNuc','Periphery','Nuclear Interior'])
ctrl_df['Treatment']='Untreated'
treat_df = pd.DataFrame(tR,columns=['TotalNuc','Periphery','Nuclear Interior'])
treat_df['Treatment']='Auxin'

# get data into longview format to plot the barplot all at once
totalDf = pd.concat([ctrl_df[ctrl_df.columns[1:]],treat_df[treat_df.columns[1:]]])
melted_df = totalDf.melt(id_vars='Treatment', var_name='Nuclear Region', value_name='Normalized Intensity')

# use sns barplot 
colors = [(1, 1, 0), (1,0,1)]*(len(melted_df)//2)
sns.set_style('whitegrid')
matplotlib.rcParams["figure.dpi"] =  100
fig,axs = plt.subplots(1,1)

sns.barplot(x='Nuclear Region', y='Normalized Intensity', alpha=0.75,hue='Treatment'
           ,data=melted_df,ax=axs,palette=colors,capsize=0.025)

axs.set_ylim(0,1.5) # edit this when your values are greater than 1.5 
axs.tick_params(axis='both', which='major', labelsize=16)
axs.set_xlabel('Nuclear Region',fontsize=20),axs.set_ylabel('Normalized Intensity',fontsize=20)
axs.legend(fontsize=16,bbox_to_anchor=(1.05,1.05),ncol=1,fancybox=True)
#plt.savefig(os.path.join(figuresPath,f'NSI_Results_Barplots_iter{5}_v3.png'),transparent=True,bbox_inches='tight')
plt.show() 


# t-test 
# get only the categories 
from scipy.stats import ttest_ind

p_df = melted_df[melted_df['Nuclear Region']=='Periphery']
i_df = melted_df[melted_df['Nuclear Region']=='Interior']

# do ttest
periphery_ttest = ttest_ind(p_df[p_df['Treatment']=='Auxin']['Normalized Intensity'],
                            p_df[p_df['Treatment']=='Untreated']['Normalized Intensity'])

interior_ttest = ttest_ind(i_df[i_df['Treatment']=='Auxin']['Normalized Intensity'],
                           i_df[i_df['Treatment']=='Untreated']["Normalized Intensity"])

print(periphery_ttest,interior_ttest) # no significnace 