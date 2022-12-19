import os
from glob import glob
from pickle import FALSE
import pandas as pd 
import cv2 
import numpy as np
from joblib import Parallel, delayed 
import joblib
import sys

def chromo_cut_off(cell_data,frames):
    area = [0] * len(frames)
    chromo_cut = [0] * len(frames)
    for i in range(len(frames)):
        temp = cell_data[cell_data.frame == frames[i]]
        f_area = 0 
        for j in temp['id'].unique():
            temp1 = temp[temp.id == j]
            x,y = temp1["x"],temp1["y"]
            cnt = contourBuilder(x,y)
            f_area = f_area + int(cv2.contourArea(cnt))
        area[i] = f_area
        chromo_cut[i] = frames[i]
    if len(area) == 0:
        return 0
    max_chromo_size = np.amax(area)
    idx = np.where(area == max_chromo_size)[0][0]
    chromo_cut = np.delete(chromo_cut, list(range(0,idx)))
    area = np.delete(area, list(range(0,idx)))
    for i in range(len(area)):
        if ((area[i] < (max_chromo_size*0.1)) & (area[i] > 250)): #note the extra addition helps to avoid gaps in the data
            return chromo_cut[i]
    return max(frames)

def contourBuilder(x,y):
    size = len(x)
    results = np.zeros(shape=(size,1,2)).astype(np.int32)
    idx = range(0,size)
    for i,j,k in zip(idx,x,y):
        results[i]=[j,k]
    return results

def max_area_finder(cell_data,frames):
    area = 0
    exptrapolation_frame = 0
    for i in frames:
        temp = cell_data[cell_data.frame == i]
        x,y = temp["x"],temp["y"]
        cnt = contourBuilder(x,y)
        temp_area = cv2.contourArea(cnt)
        if temp_area > area:
            area = temp_area
            exptrapolation_frame = i
    return exptrapolation_frame

def unpacking_data(t,xy_pos):
    contour_dataframe = np.empty((0,6))
    for i in t:
        contour_dataframe = np.vstack([contour_dataframe,i])
    xy_l = np.full((contour_dataframe.shape[0],1),xy_pos)
    contour_dataframe = np.concatenate((contour_dataframe,xy_l),1)
    df = pd.DataFrame(contour_dataframe, columns = ['x','y','id','channel','cell_id','frame',"xy"])
    return df

def fl_contour_match(phase_cnt,fl_frame):
    fl_im = cv2.imread(fl_frame)
    fl_im = np.amax(fl_im, axis=2)
    fl_im[fl_im > 0]=255
    mask = np.zeros((2424,2424),np.uint8)
    c_img = cv2.drawContours(mask, [phase_cnt], -1, (255),-1)
    fl_loc = np.where(fl_im == c_img, fl_im,0)
    kernel = np.ones((3,3),np.uint8)
    fl_loc = cv2.dilate(fl_loc, kernel, iterations = 1)
    nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(fl_loc, connectivity=8)
    fl_count = 1
    fl_info = np.empty((0,3))
    for i in range(1,nb_components):
        cnt = contour_calc(output,i,0)
        x, y = cnt[0][:,0,0], cnt[0][:,0,1]
        id = [fl_count] * len(cnt[0])
        temp = np.column_stack([x,y,id])
        fl_info = np.vstack([fl_info,temp])
        fl_info = fl_info.astype(int)
        fl_count = fl_count + 1 
    return fl_info

def find_max_list_idx(list):
    list_len = [len(i) for i in list]
    return np.argmax(np.array(list_len))

def contour_calc(mask,id,phase):
    blank = np.zeros((mask.shape[0],mask.shape[1]))
    blank[mask == (id + phase)] = 255
    blank = blank.astype(np.uint8)
    kernel = np.ones((3,3),np.uint8)
    dilated_img = cv2.dilate(blank, kernel, iterations = 2)
    contour, hierachy = cv2.findContours(dilated_img, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    return contour

def dist(gtx,gty,px,py):
    distance = (((gtx - px)**2 + (gty-py)**2)**0.5)
    min_value = min(distance)
    return np.where(distance == min_value)[0][0]

def phase_contour(cell_x,cell_y,col_mask):
    mask_val = np.unique(col_mask)
    min_mask_val = np.amin(mask_val)
    col_mask[col_mask > min_mask_val]=255
    col_mask[col_mask <= min_mask_val]=0
    nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(col_mask, connectivity=8)
    cnt_id = dist(cell_x,cell_y,centroids[1:,0],centroids[1:,1])
    cnt = contour_calc(output, cnt_id,1)
    kx,ky = cnt[0][:,0,0], cnt[0][:,0,1]
    return cnt[0], kx, ky

def single_cell_phase(x,y,grey_phase,fl,xy_pos,idx):
    cell_cnt, x_cnt,y_cnt = phase_contour(x,y,grey_phase)
    sc_info = np.column_stack([x_cnt,y_cnt,[1]*len(x_cnt),[phase_channel]*len(x_cnt)])
    if len(fl) != 0:
        for j in fl:
            fl_frame = "./"+str(j)+"/xy"+str(xy_pos)+"/raw_mask_avg/"+str(idx)+".png"     
            fl_cnt = fl_contour_match(cell_cnt,fl_frame)
            fl_chan = np.full((fl_cnt.shape[0],1),j)
            fl_holder = np.concatenate((fl_cnt,fl_chan),1)
            sc_info = np.vstack([sc_info,fl_holder])
    return sc_info
            
def cell_finder_phase(xy_pos,idx,phase_frame,t_data,fl):
    df = t_data[t_data.frame == idx]
    df = df.reset_index()
    phase_im = cv2.imread(phase_frame)
    gray = np.amax(phase_im, axis=2)
    frame_info = np.empty((0,6))
    for i in range(len(df)):
        x = df["x"][i]
        y = df["y"][i]
        sc = single_cell_phase(x,y,gray,fl,xy_pos,idx)
        cell = np.full((sc.shape[0],1),int(df["particle"][i]))
        frame = np.full((sc.shape[0],1),idx)
        d = np.concatenate((sc,cell,frame),1)
        frame_info = np.vstack([frame_info,d])
    return frame_info
        
def extrap_processing(df,image_length,p_list):
    cell_id_list = df['cell_id'].unique()
    extrapolated_contour_dataframe = np.empty((0,3))
    for i in range(len(cell_id_list)):
        if (i in p_list) == False: 
            f = df[df.cell_id == cell_id_list[i]]["frame"].unique()
            max_f = np.max(f)
            ref_frame = max_area_finder(df[df.cell_id == cell_id_list[i]],f)
            extrap_frames = list(range(max_f+1, image_length))
            if len(extrap_frames) != 0:
                extrap_cell_id = np.full((len(extrap_frames),1),cell_id_list[i])
                extrap_ref_frame = np.full((len(extrap_frames),1),ref_frame)
                a = np.vstack(extrap_frames)
                d = np.concatenate((np.vstack(extrap_frames),extrap_cell_id,extrap_ref_frame),1)
                extrapolated_contour_dataframe = np.vstack([extrapolated_contour_dataframe,d])
    return extrapolated_contour_dataframe.astype(int)

def extrapolate_checker(xy_pos,extrap_pos,df,fl):
    target_frame = int(df[extrap_pos,0])
    c_id = df[extrap_pos,1]
    ref_frame = df[extrap_pos,2]
    cnt_x, cnt_y = xy_data_phase[(xy_data_phase.cell_id == c_id) & (xy_data_phase.frame == ref_frame)]["x"].tolist(),xy_data_phase[(xy_data_phase.cell_id == c_id) & (xy_data_phase.frame == ref_frame)]["y"].tolist()
    phase_cnt = contourBuilder(cnt_x,cnt_y)
    sc_info = np.empty((0,6))
    for i in fl:
        fl_frame = "./"+str(i)+"/xy"+str(xy_pos)+"/raw_mask_avg/"+str(target_frame)+".png"
        fl_cnt = fl_contour_match(phase_cnt,fl_frame)
        fl_chan = np.full((fl_cnt.shape[0],1),i)
        fl_id = np.full((fl_cnt.shape[0],1),c_id)
        fl_frame = np.full((fl_cnt.shape[0],1),target_frame)
        fl_holder = np.concatenate((fl_cnt,fl_chan,fl_id,fl_frame),1)
        sc_info = np.vstack([sc_info,fl_holder])
    return sc_info


####################################################################
#parameters that need to be set
os.chdir(sys.argv[1])
phase_channel = 'c3'
chromosome_channel = "c2"
protein_channel = "c1"
#xy = range(1,2)#note that this can be automated for batch processing 
xy = [1]#note that this can be automated for batch processing 
####################################################################

channels = set(sorted(glob("*"))) - set(glob("*.txt"))
fl_channels = channels
fl_channels.remove(phase_channel)
fl_channels.remove(protein_channel)

for i in xy:
    tracked_data = pd.read_csv("./"+str(phase_channel)+"/xy"+str(i)+"/unlabeled/Masks_per_frame/tracks_update.csv")
    tracked_data.apply(pd.to_numeric)
    tracked_data.sort_values(["frame","particle"],ascending=[True,True])
    pngs = sorted(glob("./"+str(phase_channel)+"/xy"+str(i)+"/unlabeled/Masks_per_frame/*.png"))
    contour_dataframe = np.empty((0,6))
    num_cores = joblib.cpu_count()
    t = Parallel(n_jobs = num_cores - 2, backend="threading",verbose=5)(delayed(cell_finder_phase)(i,index,j,tracked_data,fl_channels) for index, j in enumerate(pngs))
    xy_data = unpacking_data(t,i)
    cols=[k for k in xy_data.columns if k not in ["channel"]]
    for col in cols:
        xy_data[col]=pd.to_numeric(xy_data[col])
    xy_data_phase = xy_data[xy_data.channel == phase_channel]
    if len(fl_channels) != 0:
        parental_ids = tracked_data['parent'].unique()
        extrap_data = extrap_processing(xy_data_phase,len(pngs), parental_ids)
        t = Parallel(n_jobs = num_cores - 2, backend="threading",verbose=5)(delayed(extrapolate_checker)(i,j,extrap_data,fl_channels) for j in range(len(extrap_data)))
        extrap_xy_data = unpacking_data(t,i)
        
    df = pd.concat([xy_data, extrap_xy_data])
    df.sort_values(by=['frame', 'cell_id','channel','id'])
    cols=[k for k in df.columns if k not in ["channel"]]
    for col in cols:
        df[col]=pd.to_numeric(df[col])
    df["parent"] = 0
    df["ancestor"] = 0
    for j in df["cell_id"].unique():
        if (j in parental_ids) == False:
            df.loc[df["cell_id"] == j, "parent"] = tracked_data[tracked_data.particle == j]["parent"].tolist()[0]
            temp = df[(df.cell_id == j) & (df.channel == "c2")]
            chromo_final_frame = chromo_cut_off(temp,temp["frame"].unique())
            if chromo_final_frame != 0:
                if np.max(temp['frame'].unique()) > chromo_final_frame:
                    df = df.drop(df[((df.cell_id == j) & (df.frame > chromo_final_frame))].index)
                    df = df.reset_index(drop=True)
        else:
            df.loc[df["cell_id"] == j, "ancestor"] = 1 #this means that the filament HAS NOT divided --- thus it IS an ancestor cell 
    df.to_csv("./"+str(phase_channel)+"/xy"+str(i)+"/contour_output_xy"+str(i)+".csv",index=False)
