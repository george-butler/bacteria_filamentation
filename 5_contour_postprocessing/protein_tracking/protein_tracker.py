import os
from glob import glob
from pickle import TRUE
import tracemalloc
import pandas as pd 
import cv2 
import numpy as np
from joblib import Parallel, delayed 
import joblib
import sys

def contourBuilder(x,y):
    size = len(x)
    results = np.zeros(shape=(size,1,2)).astype(np.int32)
    idx = range(0,size)
    for i,j,k in zip(idx,x,y):
        results[i]=[j,k]
    return results

def unpacking_data(t,xy_pos):
    contour_dataframe = np.empty((0,7))
    for i in t:
        contour_dataframe = np.vstack([contour_dataframe,i])
    xy_l = np.full((contour_dataframe.shape[0],1),xy_pos)
    contour_dataframe = np.concatenate((contour_dataframe,xy_l),1)
    df = pd.DataFrame(contour_dataframe, columns = ['x','y','id','channel','cell_id','frame',"tracked_protein","xy"])
    return df

def contour_calc(mask,id,correction):
    blank = np.zeros((mask.shape[0],mask.shape[1]))
    blank[mask == (id + correction)] = 255
    blank = blank.astype(np.uint8)
    contour, hierachy = cv2.findContours(blank, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    return contour

def dist(gtx,gty,px,py):
    distance = (((gtx - px)**2 + (gty-py)**2)**0.5)
    min_value = min(distance)
    return np.where(distance == min_value)[0][0]

def contour(cell_x,cell_y,col_mask):
    mask_val = np.unique(col_mask)
    min_mask_val = np.amin(mask_val)
    col_mask[col_mask > min_mask_val]=255
    col_mask[col_mask <= min_mask_val]=0
    nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(col_mask, connectivity=8)
    cnt_id = dist(cell_x,cell_y,centroids[1:,0],centroids[1:,1])
    cnt = contour_calc(output, cnt_id,1)
    kx,ky = cnt[0][:,0,0], cnt[0][:,0,1]
    return cnt[0], kx, ky

def single_cell(x,y,grey_img):
    cell_cnt, x_cnt,y_cnt = contour(x,y,grey_img)
    sc_info = np.column_stack([x_cnt,y_cnt,[1]*len(x_cnt),[protein_channel]*len(x_cnt)])
    return sc_info
            
def cell_finder(idx,target_frame,t_data):
    df = t_data[t_data.frame == idx]
    df = df.reset_index()
    im = cv2.imread(target_frame)
    gray = np.amax(im, axis=2)
    frame_info = np.empty((0,7))
    for i in range(len(df)):
        x = df["x"][i]
        y = df["y"][i]
        sc = single_cell(x,y,gray)
        cell = np.full((sc.shape[0],1),int(df["phase_id"][i]))
        frame = np.full((sc.shape[0],1),idx)
        tracked_protein = np.full((sc.shape[0],1),int(df["particle"][i]))
        d = np.concatenate((sc,cell,frame,tracked_protein),1)
        frame_info = np.vstack([frame_info,d])
    return frame_info

def contour_matcher():
    global tracked_data
    for val in range(len(tracked_data)):
        t_x,t_y = tracked_data["x"][val],tracked_data["y"][val]
        f_data = phase_contour_data[phase_contour_data.frame == tracked_data.frame[val]]
        for i in f_data['cell_id'].unique():
            x_phase,y_phase = f_data.loc[f_data["cell_id"] == i, "x"], f_data.loc[f_data['cell_id'] == i, "y"]
            cnt = contourBuilder(x_phase,y_phase)
            d = cv2.pointPolygonTest(cnt,(t_x,t_y),False)
            if d >= 0:
                tracked_data.loc[val,"phase_id"] = i

def most_frequent(List):
    counter = 0
    num = List[0]
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
    return num
    
def protein_trimmer():
    global tracked_data
    for i in tracked_data["particle"].unique():
        lst = tracked_data.loc[tracked_data["particle"] == i, "phase_id"].tolist()
        bg_check = all(elment == 0 for elment in lst) 
        if bg_check == True:
            tracked_data = tracked_data.drop(tracked_data[tracked_data.particle == i].index)
            tracked_data = tracked_data.reset_index(drop=True)
        else:
            lst = [val for val in lst if val != 0]
            m_fqt = most_frequent(lst)
            for j in np.unique(lst):
                if j != m_fqt:
                    tracked_data = tracked_data.drop(tracked_data[(tracked_data.particle == i) & (tracked_data.phase_id == j)].index)
                    tracked_data = tracked_data.reset_index(drop=True)
####################################################################
#parameters that need to be set
os.chdir(sys.argv[1])
protein_channel = "c1"
chromosome_channel = "c2"
phase_channel = 'c3'
#xy = range(1,2)#note that this can be automated for batch processing 
xy = [1]
####################################################################

for i in xy:
    contour_data = pd.read_csv("./"+str(phase_channel)+"/xy"+str(i)+"/contour_output_xy"+str(i)+".csv")
    cols=[k for k in contour_data.columns if k not in ["channel"]]
    for col in cols:
        contour_data[col]=pd.to_numeric(contour_data[col])
    phase_contour_data = contour_data[contour_data.channel == phase_channel]
    tracked_data = pd.read_csv("./"+str(protein_channel)+"/xy"+str(i)+"/unlabeled/tracks.csv")
    tracked_data.apply(pd.to_numeric)
    tracked_data.sort_values(["frame","particle"],ascending=[True,True])
    tracked_data["phase_id"] = 0
    contour_matcher()
    protein_trimmer()
    pngs = sorted(glob("./"+str(protein_channel)+"/xy"+str(i)+"/unlabeled/Masks_per_frame/*.png"))
    contour_dataframe = np.empty((0,7))
    num_cores = joblib.cpu_count()
    t = Parallel(n_jobs = num_cores - 2, backend="threading",verbose=5)(delayed(cell_finder)(index,j,tracked_data) for index, j in enumerate(pngs))
    xy_data = unpacking_data(t,i)
    cols=[k for k in xy_data.columns if k not in ["channel"]]
    for col in cols:
        xy_data[col]=pd.to_numeric(xy_data[col])
    #check if any proteins have been included that are never "in" a cell 
    for j in xy_data["tracked_protein"].unique():
        c_id_list = xy_data.loc[xy_data["tracked_protein"] == j, "cell_id"].tolist()
        if (len(c_id_list) == 1) & (c_id_list[0] == 0):
            xy_data = xy_data.drop(xy_data[xy_data["tracked_protein"] == j].index)
            xy_data = xy_data.reset_index(drop = True)
    xy_data["parent"] = 0
    xy_data["ancestor"] = 0
    for j in phase_contour_data["cell_id"].unique():
        max_f = np.amax(phase_contour_data.loc[phase_contour_data["cell_id"] == j, "frame"])
        min_f = np.amin(phase_contour_data.loc[phase_contour_data["cell_id"] == j, "frame"])
        for k in np.unique(xy_data.loc[xy_data["cell_id"] == j, "tracked_protein"].tolist()):
            xy_data = xy_data.drop(xy_data[((xy_data.tracked_protein == k) & (xy_data.frame > max_f))].index)
            xy_data = xy_data.drop(xy_data[((xy_data.tracked_protein == k) & (xy_data.frame < min_f))].index)
            xy_data = xy_data.reset_index(drop=True)
        cell_parent = phase_contour_data.loc[((phase_contour_data["cell_id"] == j) & (phase_contour_data["frame"] == max_f)), "parent"].tolist()
        cell_ancestor = phase_contour_data.loc[((phase_contour_data["cell_id"] == j) & (phase_contour_data["frame"] == max_f)), "ancestor"].tolist()
        xy_data.loc[xy_data["cell_id"] == j, "parent"] = cell_parent[0]
        xy_data.loc[xy_data["cell_id"] == j, "ancestor"] = cell_ancestor[0]

    contour_data["tracked_protein"] = 0
    df = pd.concat([xy_data, contour_data], sort=False)
    df.sort_values(by=['frame', 'cell_id','channel','id'])
    cols=[k for k in df.columns if k not in ["channel"]]
    for col in cols:
        df[col]=pd.to_numeric(df[col])
    df.to_csv("./"+str(phase_channel)+"/xy"+str(i)+"/contour_output_xy"+str(i)+"_with_protein_tracking.csv",index=False)

#x and y are the coordinates of the contour
