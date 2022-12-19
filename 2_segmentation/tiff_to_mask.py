from cmath import phase
from joblib import Parallel, delayed 
import joblib
import cv2
import numpy as np
import os
from glob import glob 
import sys 

def remove_small_objects(img, min_size, max_size):
    nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(img, connectivity=8)
    sizes = stats[1:, -1]
    nb_components = nb_components - 1
    img2 = img
    for i in range(0, nb_components):
        if (sizes[i] < min_size) | (sizes[i] > max_size):
            img2[output == i + 1] = 0

    return img2 

def invert_check(fig):
    blur = cv2.GaussianBlur(fig,(3,3),0)
    thresh = cv2.adaptiveThreshold(blur, 255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,97,2)
    thresh = cv2.bitwise_not(thresh)
    no_black_pix = np.sum(thresh == 0)
    invert = 255 - fig
    blur_inv = cv2.GaussianBlur(invert,(3,3),0)
    thresh_inv = cv2.adaptiveThreshold(blur_inv, 255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,97,2)
    thresh_inv = cv2.bitwise_not(thresh_inv)
    no_black_pix_inv = np.sum(thresh_inv == 0)
    if no_black_pix_inv > no_black_pix:
        return thresh_inv
        
    else:
        return thresh

def seg_machine(fig,channel):
    if channel == protein_channel:
        ret2,thresh = cv2.threshold(fig,120,255,cv2.THRESH_BINARY)
        thresh = remove_small_objects(thresh,20,1000)
        ret, markers = cv2.connectedComponents(thresh)
        markers = markers.astype(np.uint8)
        return markers
    if channel == chromosome_channel:
        ret2,thresh = cv2.threshold(fig,80,255,cv2.THRESH_BINARY)
        ret, markers = cv2.connectedComponents(thresh)
        markers = markers.astype(np.uint8)
        return markers
    else:
        thresh = invert_check(fig)
        mask1 = np.zeros((thresh.shape[0]+2, thresh.shape[0]+2), np.uint8) 
        output = cv2.floodFill(thresh,mask1,(0,0),255)
        mask_inv=cv2.bitwise_not(thresh)
        kernel = np.ones((3,3),np.uint8)
        opening = cv2.morphologyEx(mask_inv,cv2.MORPH_OPEN,kernel, iterations = 1)
        sure_bg = cv2.dilate(opening,kernel,iterations=5)
        dist_transform = cv2.distanceTransform(opening,cv2.DIST_L2,5)
        ret, sure_fg = cv2.threshold(dist_transform,0.05*dist_transform.max(),255,0)
        sure_fg = np.uint8(sure_fg)
        unknown = cv2.subtract(sure_bg,sure_fg)
        sure_fg = remove_small_objects(sure_fg,(thresh.shape[0] / 10), ((thresh.shape[0]**2)/2))
        ret, markers = cv2.connectedComponents(sure_fg)
        markers[unknown==255] = 0
        markers = markers.astype(np.uint8)
    return markers


def file_management(idx,xy_path,channel):
        img = cv2.imread(xy_path + "/raw/" + str(idx)+".tif")
        seg_output = seg_machine(cv2.cvtColor(img,cv2.COLOR_BGR2GRAY),channel)
        cv2.imwrite(xy_path + "/raw_mask_avg/" + str(idx) + ".png", seg_output)

#################################################################################
path = sys.argv[1]
phase_channel = 'c3'
chromosome_channel = 'c2'
protein_channel = 'c1'
os.chdir(path)

seperated_files = os.path.isdir(path + "/images")
if seperated_files == "False":
    print("No seperated .tifs")
    sys.exit()

phase_channel_folder = os.path.isdir(path + "/images/" + phase_channel)
if phase_channel_folder == "False":
    print("No phase channel")
    sys.exit()

os.chdir(path + "/images/")
channels = set(sorted(glob("*"))) - set(glob("*.csv"))
channels.remove(phase_channel)

xy_phase = glob("./" + phase_channel + "/*")
for i in xy_phase:
    os.mkdir(i + "/raw_mask_avg")
    no_phase_img = len(glob(i + "/raw/*.tif"))
    num_cores = joblib.cpu_count()
    Parallel(n_jobs= num_cores - 2, backend = "threading", verbose=5)(delayed(file_management)(j,i,phase_channel) for j in range(no_phase_img))

if len(channels) != 0:
    for i in channels:
        xy_fl = glob("./" + i + "/*")
        for j in xy_fl:
            os.mkdir(j + "/raw_mask_avg")
            no_fl_img = len(glob(j + "/raw/*.tif"))
            num_cores = joblib.cpu_count()
            Parallel(n_jobs= num_cores - 2, backend = "threading", verbose=5)(delayed(file_management)(k,j,i) for k in range(no_fl_img))

