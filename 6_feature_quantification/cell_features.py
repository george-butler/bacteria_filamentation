import pandas as pd 
import numpy as np 
import cv2 
from joblib import Parallel, delayed 
import joblib
from skimage.morphology import medial_axis
import sys

def unpacking_data(t,xy_pos):
    contour_dataframe = np.empty((0,t[0].shape[1]))
    for i in t:
        contour_dataframe = np.vstack([contour_dataframe,i])
    xy_l = np.full((contour_dataframe.shape[0],1),xy_pos)
    contour_dataframe = np.concatenate((contour_dataframe,xy_l),1)
    df = pd.DataFrame(contour_dataframe, columns = ['frame','cell_id','channel','id','area','perimeter','medial_axis',"parent","ancestor","xy"])
    return df

def contourBuilder(x,y):
	size = len(x)
	results = np.zeros(shape = (size,1,2)).astype(np.int32)
	idx = range(0,size)
	for i,j,k in zip(idx,x,y):
		results[i] = [j,k]
	return results

def med_axis_calc(cnt):
    mask = np.zeros((2424,2424),np.uint8)
    img = cv2.drawContours(mask, [cnt], -1, (255),-1)
    img[img > 0]=255
    skel = medial_axis(img)
    return len(skel[skel == True])


def measure(x,y):
    cell_contour = contourBuilder(x,y)
    area = int(cv2.contourArea(cell_contour))
    perimeter = int(cv2.arcLength(cell_contour, True))
    med_axis = med_axis_calc(cell_contour)
    return area, perimeter, med_axis

def frame_analysis(f,d):
    frame_info = np.empty((0,9))
    frame_data = d[d.frame == f]
    for i in frame_data["cell_id"].unique():
        temp = frame_data[frame_data.cell_id == i]
        for j in temp["channel"].unique():
            temp1 = temp[temp.channel == j]
            for k in temp1["id"].unique():
                a,p,ma = measure(temp1[temp1.id == k]['x'],temp1[temp1.id == k]['y'])
                parent = temp1[temp1.id == k]["parent"].tolist()[0]
                ancestor = temp1[temp1.id == k]["ancestor"].tolist()[0]
                d = [f,i,j,k,a,p,ma,parent,ancestor]
                frame_info = np.vstack([frame_info,d])
    return frame_info

xy = [1]
protein_channel = "c1"
for k in xy:		
    data = pd.read_csv(sys.argv[1] + "/c3/xy"+str(k)+"/" + "contour_output_xy"+str(k)+"_with_protein_tracking.csv")
    data = data[data.channel != protein_channel]
    num_cores = joblib.cpu_count()
    f_list = data["frame"].unique()
    t = Parallel(n_jobs = num_cores - 14, backend="threading",verbose=5)(delayed(frame_analysis)(i,data) for i in f_list)

    xy_data = unpacking_data(t,k)
    xy_data.to_csv(sys.argv[1] + "/c3/xy"+str(k)+"/" + "cell_features_xy"+str(k)+".csv", index = False)
