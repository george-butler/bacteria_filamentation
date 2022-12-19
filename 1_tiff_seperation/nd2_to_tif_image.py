import itertools
from joblib import Parallel, delayed 
import joblib
from nd2file import ND2MultiDim
import cv2
import numpy as np
import os
import glob 
import sys 
import datetime

def ripper(file_name, index):
    img16 = file_name.image_singlechannel(multipoint=index[2],timepoint=index[1],channel=index[0])
    ratio = np.amax(img16) / 255
    img8 = (img16 / ratio).astype('uint8')
    cv2.imwrite(directory_location + "/c" + str(index[0]+1) + "/xy" + str(index[2]+1)+ "/" + str(index[1]) + ".tif", img8)

path = sys.argv[1]

os.chdir(path)

original_files = glob.glob("*.nd2")

if len(original_files) == 0:
    print("No files in this directory")
    sys.exit()

for file in glob.glob("*.nd2"):
    if os.path.exists(path+"/images") == True:
        print("Warning this file may have already been split")
        directory_location = path + "/images" + datetime.datetime.now().isoformat()
        os.mkdir(directory_location)
    else:
        directory_location = path + "/images"
        os.mkdir(directory_location)

    nd2 = ND2MultiDim(path + "/" + file)

    no_xy = nd2.multipointcount
    no_chan = nd2.channels
    no_timepoints = nd2.timepointcount


    lines = ['Metadata', "Calibration " + str(nd2.calibration),
    "Number of points " + str(no_xy), "Number of channels " + str(no_chan), "Number of timepoints " + str(no_timepoints),
    "Image Height " + str(nd2.height), "Image Width " + str(nd2.width)]
    with open(path + "/images" + "/image_aquisition_information.txt", 'w') as f:
        f.write('\n'.join(lines))

    lxy = list(range(no_xy))
    lc = list(range(no_chan))
    for i in lc:
        os.mkdir(directory_location + "/c" + str(i+1))
        for j in lxy:
            os.makedirs(directory_location + "/c" + str(i+1) + "/xy" + str(j+1) + "/raw")
            
    ltimepoint = list(range(no_timepoints))
    list_combo = list(itertools.product(lc,ltimepoint,lxy))
    num_cores = joblib.cpu_count()
    print("Total number of tasks to complete",len(list_combo))
    Parallel(n_jobs= num_cores - 12, backend = "threading", verbose=5)(delayed(ripper)(nd2,i) for i in list_combo)
