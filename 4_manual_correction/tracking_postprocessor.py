import pandas as pd 
import numpy as np 
import collections
import sys 
#this function checks whether the cell has another link further on in time 
#it also ensure that any change in id doesn't include the current frame 
def cell_checker(dataframe, value,set_value):
    new_id_list = df.iloc[:,3].tolist()
    comparison_value = df.iloc[value][1]
    if (comparison_value in new_id_list):
        h = [i for i,val in enumerate(new_id_list) if val == df.iloc[value][1]]
        idx = None
        for i in h:
            if df.iloc[i][0] > df.iloc[value][0]:
                idx = i
                break
        return idx
    else:
        return None
#changes values in the data frame 
def cell_changer(dataframe,idx, set_value):
    df.iloc[idx,4] = set_value
    df.iloc[idx,3] = -1
    return idx
#works through the dataframe and checks if frame has already been considered 
def main(data):
    df["corrected_id"] = 0
    for i in range(len(df)):
        if df.iloc[i][3] != -1:
            set_value = df.iloc[i][3]
            df.iloc[i,4] = set_value
            df.iloc[i,3] = -1
            output = True
            value = i
            while output == True:
                index = cell_checker(df, value,set_value)
                if index != None:
                    value = cell_changer(df,index, set_value)
                else:
                    output = False
    return df.sort_values(["corrected_id","frame_id"], ascending = [True,True])

t_path = sys.argv[1]

tracks = pd.read_csv(t_path + str("/unlabeled/tracks.csv"))
tracks.apply(pd.to_numeric)

notes_data = pd.read_csv(t_path + str("/manual_notes.csv"))
notes_data.apply(pd.to_numeric)

notes_data["frame_id"] = notes_data["frame_id"] - 1

lineage_date = notes_data[(notes_data.key == 1)]
lineage_date = lineage_date.reset_index(drop = True)
lineage_date["updated_id_parent"] = 0
div_number = len(lineage_date)

a = notes_data["key"].unique()

#1 records lineage information 
#2 corrects cases where cells labels have been mixed up
#3 delete rows that contain a given cell id and THE frame id 
#4 delete rows that contain a give cell id and ALL FURTHER frame ids
#5 delete all rows that contain a given cell id 

#deletes all rows that contain a given cell id 
if 5 in a:
    all_row_del = notes_data[notes_data.key == 5]
    all_row_del = all_row_del.reset_index(drop = True)
    for index, row in all_row_del.iterrows():
        tracks = tracks.drop(tracks[tracks.particle == row["cell_id"]].index)
    tracks = tracks.reset_index(drop = True)
#deletes the rows that contain a given cell id and the frame id 
if 4 in a:
    row_del = notes_data[notes_data.key == 4]
    row_del = row_del.reset_index(drop = True)
    for index, row in row_del.iterrows():
        tracks = tracks.drop(tracks[(tracks.particle == row["cell_id"]) & (tracks.frame >= row["frame_id"])].index)
    tracks = tracks.reset_index(drop = True)

if 3 in a:
    row_del = notes_data[notes_data.key == 3]
    row_del = row_del.reset_index(drop = True)
    for index, row in row_del.iterrows():
        tracks = tracks.drop(tracks[(tracks.particle == row["cell_id"]) & (tracks.frame == row["frame_id"])].index)
    tracks = tracks.reset_index(drop = True)


if 2 in a:
    df = notes_data[notes_data.key == 2]
    df = df.reset_index(drop = True)
    corrected_label_data = main(df)
    corrected_id_uni = corrected_label_data["corrected_id"].unique()
    corrected_label_data = corrected_label_data.sort_values("frame_id", ascending = True)
    tracks["corrected_id"] = tracks["particle"]
    for i in corrected_id_uni:
        temp = corrected_label_data[corrected_label_data.corrected_id == i]
        temp = temp.reset_index(drop = True)
        if len(temp) == 1:
            for i in range(len(temp)):
                tracks.loc[(tracks.particle == temp["cell_id"][i]) & (tracks.frame >= temp["frame_id"][i]), "corrected_id"] = temp["corrected_id"][i]
                if div_number != 0:
                    lineage_date.loc[(lineage_date.cell_id == temp["cell_id"][i]) & (lineage_date.frame_id >= temp["frame_id"][i]), "updated_id_parent"] = temp["corrected_id"][i]
        else:
            for j in range(len(temp)):
                if j == (len(temp) - 1):
                    tracks.loc[(tracks.particle == temp["cell_id"][j]) & (tracks.frame >= temp["frame_id"][j]), "corrected_id"] = temp["corrected_id"][j]
                    if div_number != 0:
                        lineage_date.loc[(lineage_date.cell_id == temp["cell_id"][j]) & (lineage_date.frame_id >= temp["frame_id"][j]), "updated_id_parent"] = temp["corrected_id"][j]
                else:
                    tracks.loc[(tracks.particle == temp["cell_id"][j]) & (tracks.frame >= temp["frame_id"][j]) & (tracks.frame < temp["frame_id"][(j+1)]), "corrected_id"] = temp["corrected_id"][j]
                    if div_number != 0:
                        lineage_date.loc[(lineage_date.cell_id == temp["cell_id"][j]) & (lineage_date.frame_id >= temp["frame_id"][j]) & (lineage_date.frame_id < temp["frame_id"][(j+1)]), "updated_id_parent"] = temp["corrected_id"][j]
    tracks["particle"] = tracks["corrected_id"]
    tracks.drop('corrected_id', axis=1, inplace=True)
    tracks = tracks.sort_values(["particle","frame"],ascending = [True,True]) 
    tracks = tracks.reset_index(drop = True)

lineage_date = lineage_date.reset_index(drop = True)

lineage_date["screen_id"] = lineage_date["new_id"]
lineage_date["tracked_id"] = lineage_date["new_id"]
if 1 in a:
    max_cell_value = max(tracks["particle"].unique())
    for index, row in lineage_date.iterrows():
        lineage_date.loc[index, "offspring_id"] = max_cell_value + index + 1
    for index, row in lineage_date.iterrows():
        lineage_date.loc[(lineage_date.cell_id == row["new_id"]) & (lineage_date.frame_id > row["frame_id"]), "cell_id"] = row["offspring_id"]
        if lineage_date.loc[index, "updated_id_parent"] == 0:
            lineage_date.loc[index, "updated_id_parent"] = lineage_date.loc[index, "cell_id"]
    for index, row in lineage_date.iterrows():
        if row["cell_id"] == row["new_id"]:
            if row["cell_id"] != row["updated_id_parent"]:
                if row["cell_id"] == row["screen_id"]:
                    lineage_date.loc[index, "tracked_id"] = lineage_date.loc[index, "updated_id_parent"]
    lineage_date["updated_id_parent"] = lineage_date["cell_id"]

    for i in lineage_date["tracked_id"].unique():
        temp = lineage_date[lineage_date.tracked_id == i]
        temp = temp.reset_index(drop = True)
        if len(temp) == 1:
            for j in range(len(temp)):
                tracks.loc[(tracks.particle == temp["tracked_id"][j]) & (tracks.frame >= temp["frame_id"][j]), "particle"] = temp["offspring_id"][j]
        else:
            for j in range(len(temp)):
                if j != len(temp)-1:
                    tracks.loc[(tracks.particle == temp["tracked_id"][j]) & (tracks.frame >= temp["frame_id"][j]) & (tracks.frame < temp["frame_id"][(j+1)]), "particle"] = temp["offspring_id"][j]
                else:
                    tracks.loc[(tracks.particle == temp["tracked_id"][j]) & (tracks.frame >= temp["frame_id"][j]), "particle"] = temp["offspring_id"][j] 

tracks = tracks.sort_values(["frame","particle"],ascending = [True,True])
tracks["parent"] = 0
for index, row in lineage_date.iterrows():
    tracks.loc[tracks.particle == row["offspring_id"], "parent"] = row["updated_id_parent"]

for i in tracks["frame"].unique():
    a = tracks[tracks.frame == i].particle
    h = [item for item, count in collections.Counter(a).items() if count > 1]
    if len(h) > 0:
        for j in h:
            tracks = tracks.drop(tracks[(tracks.particle == j) & (tracks.frame == i)].index)

tracks.to_csv(t_path+str("/unlabeled/Masks_per_frame/tracks_update.csv"), index = False)
lineage_date.to_csv(t_path+str("/lineage_date.csv"), index = False)