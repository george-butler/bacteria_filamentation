import numpy as np
import pandas as pd
import trackpy as tp
import math
from collections import defaultdict, Counter
from datetime import datetime
from os.path import join, split, normpath
from skimage.measure import regionprops
from skimage.morphology import binary_dilation

tp.quiet(suppress=True)


def track_cells(path_in, raw_masks, no_apps, min_cell_id=None):
    min_number_of_app = no_apps
    n_active_features, idx_start_active_features, col_tuple, col_weights, \
        max_displacement, max_absence_interval, out_path \
        = initialize_experiment_parameters(path_in)
    trj = calculate_initial_cell_info(raw_masks, n_active_features, idx_start_active_features,
                                      col_tuple, col_weights, max_displacement, max_absence_interval)
    relink_missing_cells(raw_masks, trj, col_tuple, col_weights, idx_start_active_features, n_active_features,
                         max_absence_interval, max_displacement)
    if min_cell_id is not None:
        reindex_with_min_cell_id(trj, min_cell_id)
    a = trj['particle'].unique()
    for i in a:
        x = trj[trj.particle == i]["wtd_x"].iloc[0]
        y = trj[trj.particle == i]["wtd_y"].iloc[0]
        no_appearances = sum(trj["particle"] == i)
        if no_appearances < min_number_of_app:
            trj = trj.drop(trj[trj.particle == i].index)
        if x < 24:
            trj = trj.drop(trj[trj.particle == i].index)
        if y < 24:
            trj = trj.drop(trj[trj.particle == i].index)
        if x > 2400:
            trj = trj.drop(trj[trj.particle == i].index)
        if y > 2400:
            trj = trj.drop(trj[trj.particle == i].index)
    trj = trj.reset_index(drop=True)
    if min_cell_id is not None:
        reindex_with_min_cell_id(trj, min_cell_id)
    return trj, col_tuple, col_weights


def initialize_experiment_parameters(path_in):
    n_active_features = 8  # Must ALWAYS put the active features in a continuous!!!
    idx_start_active_features = 1  # First active feature
    col_tuple = {'original': ['frame',
                              'y', 'x', 'equivalent_diameter', 'perimeter', 'eccentricity',
                              'orientation_x_2_sin', 'orientation_x_2_cos',
                              'true_solidity',
                              'solidity',
                              'area', 'mean_intensity', 'angle', 'circularity'],
                 'extra': ['bbox_top', 'bbox_left', 'bbox_bottom', 'bbox_right']}
    # Each column will have a "sister" column with prefix 'wtd_', which will represent its weighted version
    col_tuple['weighted'] = ['wtd_{}'.format(x) for x in col_tuple['original']]
    #max_displacement = 100
    max_displacement = 25
    # Max no. of frames that a cell id may be missing (used mostly in tracking)
    max_absence_interval = 5
    col_weights = defaultdict(lambda: 1,  # Default weight is 1
                              {'frame': max_displacement / max_absence_interval,
                               'equivalent_diameter': 0.75,
                               'perimeter': 0.25,
                               'eccentricity': 5,
                               'orientation_x_2_sin': 5,
                               'orientation_x_2_cos': 5,
                               'true_solidity': 5 * np.pi})

    # Prepare output path
    in_path_part_1, in_path_part_2 = split(normpath(path_in))
    exp_time = datetime.now()  # Experiment time
    out_path = join(in_path_part_1,
                    '{}_Exp_{}-{:02d}-{:02d}T{:02d}{:02d}{:02d}'.format(in_path_part_2,
                                                                        exp_time.year,
                                                                        exp_time.month,
                                                                        exp_time.day,
                                                                        exp_time.hour,
                                                                        exp_time.minute,
                                                                        exp_time.second))
    return n_active_features, idx_start_active_features, col_tuple, col_weights, \
        max_displacement, max_absence_interval, out_path

def calculate_initial_cell_info(raw_masks, n_active_features, idx_start_active_features,
                                col_tuple, col_weights, max_displacement, max_absence_interval):
    features = pd.DataFrame()
    for i_frame in range(raw_masks.shape[2]):
        for region in regionprops(raw_masks[:, :, i_frame], intensity_image=raw_masks[:, :, i_frame]):
            if region.mean_intensity < 1:
                # Skip background (intensity 0)
                continue
            # Append all features
            features = features.append(
                [get_region_info(region, i_frame, col_weights)])
    # Compute initial tracks
    trj = tp.link_df(features,
                     search_range=max_displacement,
                     pos_columns=col_tuple['weighted'][idx_start_active_features:(
                         idx_start_active_features + n_active_features)],
                     memory=max_absence_interval,
                     neighbor_strategy='KDTree')
    # Reset indexes (current trj has index 0 for all rows; creating new indexes will be useful later)
    trj.reset_index(drop=True, inplace=True)
    return trj


def get_region_info(region, i_frame, col_weights):
    # Compute features
    feat_dict = {'y': region.centroid[0],
                 'x': region.centroid[1],
                 'equivalent_diameter': region.equivalent_diameter,
                 'perimeter': max(1, region.perimeter),
                 'eccentricity': region.eccentricity,
                 'orientation_x_2_sin': np.sin(2 * region.orientation),
                 'orientation_x_2_cos': np.cos(2 * region.orientation),
                 'true_solidity': region.equivalent_diameter / max(1, region.perimeter),
                 'solidity': region.solidity,
                 'area': region.area,
                 'mean_intensity': region.mean_intensity,
                 'angle': region.orientation,
                 'frame': i_frame,
                 'circularity' : (4 * math.pi * region.area) / (region.perimeter * region.perimeter)
                 }
    # Compute weighted features
    weighted_features_list = [('wtd_{}'.format(feat_name), col_weights[feat_name] * feat_val)
                              for feat_name, feat_val in feat_dict.items()]
    feat_dict.update(dict(weighted_features_list))
    feat_dict.update(dict([('bbox_top', region.bbox[0]),
                           ('bbox_left', region.bbox[1]),
                           ('bbox_bottom', region.bbox[2]),
                           ('bbox_right', region.bbox[3])]))
    return feat_dict


# Calculate euclidian distance between the vectors of information of the two cells
def euclid_cell_dist(cell1_info, cell2_info, col_tuple_list):
    return np.linalg.norm(
        cell1_info.loc[col_tuple_list].values - cell2_info.loc[col_tuple_list].values)


# Get index of cell info in the  trj dataframe, based on a given column value
def get_trj_idx(trj, i_frame, col_name, col_value_list):
    return trj.index[(trj['frame'] == i_frame) & (trj[col_name].isin(col_value_list))]


# Get cell neighbors and their number of pixels on the border of the given cell
def get_cell_neighbors(raw_masks, trj, i_frame, particle_id, forbidden_ids=[]):
    i_frame = int(i_frame)
    raw_mask = raw_masks[:, :, i_frame]
    cell_trj_idx = get_trj_idx(trj, i_frame, 'particle', [particle_id])
    if cell_trj_idx.size == 0:
        return []
    cell_trj_info = trj.iloc[cell_trj_idx[0]]
    # Original bounding box
    cell_bbox_orig = (cell_trj_info['bbox_top'],
                      cell_trj_info['bbox_left'],
                      cell_trj_info['bbox_bottom'],
                      cell_trj_info['bbox_right'])
    # Enlarge bounding box
    cell_bbox = (int(max(0, cell_bbox_orig[0] - 1)),
                 int(max(0, cell_bbox_orig[1] - 1)),
                 int(min(raw_mask.shape[0], cell_bbox_orig[2] + 1)),
                 int(min(raw_mask.shape[1], cell_bbox_orig[3] + 1)))
    cell_nbhood = raw_mask[cell_bbox[0]:cell_bbox[2], cell_bbox[1]:cell_bbox[3]]
    cell_mask = (
        cell_nbhood == cell_trj_info['mean_intensity']).astype(np.uint8)
    cell_border = binary_dilation(cell_mask, selem=np.ones((3, 3))) - cell_mask
    # "nb" stands for neighbor
    common_border_neighbors = Counter([x for x in cell_nbhood[cell_border > 0]
                                       if x > 0 and x not in forbidden_ids])
    # This holds all neighbor id's, ordered in descending order by length of common border
    nb_particle_ids = []
    for intensity, n_pixels in common_border_neighbors.most_common():
        nb_particle_ids.append((trj.loc[get_trj_idx(trj, i_frame, 'mean_intensity', [intensity])[0], 'particle'],
                                n_pixels))
    return nb_particle_ids


# Create new frame with the two cells merged and calculate info of the new merged cell
# prev_intensity is the mean_intensity of a potential common ancestor of p and nb in the previous frame
def merge_cell_pair(raw_masks, trj, i_frame, col_weights, p_id, nb_id):
    p_intensity = trj.loc[get_trj_idx(trj, i_frame, 'particle', [p_id])[
        0], 'mean_intensity']
    nb_intensity = trj.loc[get_trj_idx(trj, i_frame, 'particle', [nb_id])[
        0], 'mean_intensity']
    # Prepare the merged cell. All others are eliminated, in order to speed up region property calculation
    merged_frame = raw_masks[:, :, i_frame].copy()
    merged_frame[merged_frame == p_intensity] = nb_intensity
    merged_frame[merged_frame != nb_intensity] = 0
    # Calculate the region properties of the merged cell
    for region in regionprops(merged_frame, intensity_image=merged_frame):
        if region.mean_intensity == nb_intensity:
            merged_info = pd.Series(
                get_region_info(region, i_frame, col_weights))
            break
    # We create the complete frame containing the merged cells
    merged_frame = raw_masks[:, :, i_frame].copy()
    merged_frame[merged_frame == p_intensity] = nb_intensity
    return merged_frame, merged_info


# Merge cell pairs from a given frame till the end
def merge_cell_pair_sequence(raw_masks, trj, start_frame,
                             col_tuple, col_weights, cell_frames,
                             merged_frame, merged_info,
                             p_id, nb_id, target_id=None):
    n_frames = raw_masks.shape[2]
    all_cols = col_tuple['original'] + \
        col_tuple['weighted'] + col_tuple['extra']
    p_trj_idx = get_trj_idx(trj, start_frame, 'particle', [p_id])[0]
    nb_trj_idx = get_trj_idx(trj, start_frame, 'particle', [nb_id])[0]
    # For start_frame we merge anyway
    raw_masks[:, :, start_frame] = merged_frame.copy()
    trj.loc[nb_trj_idx, all_cols] = merged_info.loc[all_cols].copy()
    if target_id is not None:
        trj.loc[nb_trj_idx, 'particle'] = target_id
    # Drop row of merged cell from trajectories
    trj.drop(p_trj_idx, inplace=True)
    trj.reset_index(drop=True, inplace=True)
    # For the next frames we merge only if the two cells touch each other
    for i_frame in range(1 + start_frame, n_frames):
        flag_1 = i_frame in cell_frames[p_id]
        flag_2 = i_frame in cell_frames[nb_id]
        flag_3 = nb_id in [x[0] for x in get_cell_neighbors(
            raw_masks, trj, i_frame, p_id)]
        if (not flag_1) or (not flag_2) or (not flag_3):
            break  # The two cells do not touch each other in current frame
        p_trj_idx = get_trj_idx(trj, i_frame, 'particle', [p_id])[0]
        nb_trj_idx = get_trj_idx(trj, i_frame, 'particle', [nb_id])[0]
        # Modify nb cell
        merged_frame, merged_info = merge_cell_pair(raw_masks, trj, i_frame, col_weights,
                                                    p_id, nb_id)
        raw_masks[:, :, i_frame] = merged_frame.copy()
        trj.loc[nb_trj_idx, all_cols] = merged_info.loc[all_cols].copy()
        if target_id is not None:
            trj.loc[nb_trj_idx, 'particle'] = target_id
        # Drop p cell from trajectories
        trj.drop(p_trj_idx, inplace=True)
        trj.reset_index(drop=True, inplace=True)

# Recover cells that went missing in previous frames
# and reappeared in current frame under another particle id
def relink_missing_cells(raw_masks, trj, col_tuple, col_weights, idx_start_active_features, n_active_features,
                         max_absence_interval, max_displacement):
    for i_frame in range(1, raw_masks.shape[2]):
        relink_missing_cells_per_frame(raw_masks, trj, i_frame,
                                       col_tuple, col_weights, idx_start_active_features, n_active_features,
                                       max_absence_interval, max_displacement)


def relink_missing_cells_per_frame(raw_masks, trj, i_frame,
                                   col_tuple, col_weights, idx_start_active_features, n_active_features,
                                   max_absence_interval, max_displacement):
    all_cols = col_tuple['original'] + \
        col_tuple['weighted'] + col_tuple['extra']
    i_history_start = max(0, i_frame - max_absence_interval - 1)
    cell_frames = trj.groupby('particle')['frame'].apply(set).to_dict()
    # Get history of particles per frame
    particle_id_dict = {}
    for j_frame in range(i_history_start, 1 + i_frame):
        # We use j_frame as an alternative iterator, because i_frame comes from the caller
        particle_id_dict[j_frame] = set(
            trj[trj['frame'] == j_frame]['particle'].tolist())
    # Find newly appeared cells in current frame
    new_cell_ids = particle_id_dict[i_frame] - particle_id_dict[i_frame - 1]
    new_cell_info_dict = dict(
        [(x, trj.iloc[get_trj_idx(trj, i_frame, 'particle', [x])[0]]) for x in new_cell_ids])
    new_cell_importance_order = [id for diam, id in sorted([(cell_info['equivalent_diameter'], cell_id)
                                                            for cell_id, cell_info in new_cell_info_dict.items()],
                                                           reverse=True)]
    # Find which disappeared in the past and did not appear back
    disapp_cell_info_dict = {}
    # "Future" with respect to previous frames
    future_cell_id_list = set(particle_id_dict[i_frame])
    past_frame_idx = set(range(i_frame))
    for j_frame in range(i_frame - 1, i_history_start, -1):
        # We only consider the last disappearance of a cell during the maximum absence interval
        for disapp_cell_id in particle_id_dict[j_frame] - future_cell_id_list:
            disapp_cell_info_dict[disapp_cell_id] = trj.iloc[get_trj_idx(trj,
                                                                         j_frame,
                                                                         'particle',
                                                                         [disapp_cell_id])[0]].copy()
        future_cell_id_list |= particle_id_dict[j_frame]
    # Find best candidate for each new cell
    for new_cell_id in new_cell_importance_order:
        new_cell_info = new_cell_info_dict[new_cell_id].copy()
        candidates_needing_merge = []
        candidate_distances = []
        for disapp_cell_id, disapp_cell_info in disapp_cell_info_dict.items():
            if euclid_cell_dist(new_cell_info,
                                disapp_cell_info,
                                ['y', 'x']) > max_displacement:
                continue  # The two cells are too far away from each other
            common_frames = cell_frames[new_cell_id] & cell_frames[disapp_cell_id]
            if len(common_frames) > 0:
                if not cells_always_in_touch_in_common_frames(raw_masks, trj,
                                                              new_cell_id, disapp_cell_id, common_frames):
                    continue  # The two cells appear in the same frame either in the past or in the future
                candidates_needing_merge.append(disapp_cell_id)
            candidate_distances.append(
                (euclid_cell_dist(new_cell_info,
                                  disapp_cell_info,
                                  col_tuple['weighted'][idx_start_active_features:(idx_start_active_features
                                                                                   + n_active_features)] + ['wtd_frame']),
                 disapp_cell_id))
        if len(candidate_distances) > 0:
            # We've found something to link to
            candidate_distances.sort()
            old_cell_id = candidate_distances[0][1]
            # Change all new_cell_id into old_cell_id for trj and cell_frames
            if min(cell_frames[old_cell_id]) <= min(cell_frames[new_cell_id]):
                oldest_cell_id = old_cell_id
                newest_cell_id = new_cell_id
            else:
                oldest_cell_id = new_cell_id
                newest_cell_id = old_cell_id
            if old_cell_id in candidates_needing_merge:
                # First, we merge particles in frames where they both appear
                for common_frame in cell_frames[new_cell_id] & cell_frames[old_cell_id]:
                    common_frame = int(common_frame)
                    merged_frame, merged_info = merge_cell_pair(raw_masks, trj, common_frame, col_weights,
                                                                new_cell_id, old_cell_id)
                    raw_masks[:, :, common_frame] = merged_frame.copy()
                    old_cell_idx = get_trj_idx(
                        trj, common_frame, 'particle', [old_cell_id])[0]
                    trj.loc[old_cell_idx,
                            all_cols] = merged_info.loc[all_cols].copy()
                    # Drop new_id cell from trajectories
                    trj.drop(get_trj_idx(trj, common_frame, 'particle', [
                             new_cell_id])[0], inplace=True)
                    trj.reset_index(drop=True, inplace=True)
            # Replace newest_cell_id with oldest_cell_id in all frames
            trj.loc[trj['particle'] == newest_cell_id,
                    'particle'] = oldest_cell_id
            cell_frames[oldest_cell_id] |= cell_frames[newest_cell_id]
            del cell_frames[newest_cell_id]
            # Remove old_cell_id from being a candidate replacement for the other new cells in the same frame
            del disapp_cell_info_dict[old_cell_id]


def cells_always_in_touch_in_common_frames(raw_masks, trj, cell_id_1, cell_id_2, common_frames):
    for i_frame in common_frames:
        cell_1_info = trj.iloc[get_trj_idx(
            trj, i_frame, 'particle', [cell_id_1])[0]]
        cell_2_info = trj.iloc[get_trj_idx(
            trj, i_frame, 'particle', [cell_id_2])[0]]
        if cell_1_info['area'] < cell_2_info['area']:
            small_cell_info = cell_1_info
            large_cell_info = cell_2_info
        else:
            small_cell_info = cell_2_info
            large_cell_info = cell_1_info
        if large_cell_info['particle'] not in [x[0] for x in get_cell_neighbors(raw_masks, trj, i_frame,
                                                                                small_cell_info['particle'])]:
            return False
    return True


def reindex_with_min_cell_id(trj, min_cell_id):
    # First, move all cells to higher ids which are safe in the course of id changing
    id_offset = min_cell_id + 1 + max(trj['particle'])
    trj['particle'] = id_offset + trj['particle']
    # Now replace all ids starting from min_cell_id up
    past_cell_ids = set()
    next_cell_id = min_cell_id
    frame_cells = trj.groupby('frame')['particle'].apply(set).to_dict()
    for frame in sorted(list(set(trj['frame'].tolist()))):
        curr_cell_ids = frame_cells[frame]
        for cell_id in sorted(list(curr_cell_ids - past_cell_ids)):
            trj.loc[trj['particle'] == cell_id, 'particle'] = next_cell_id
            next_cell_id += 1
        past_cell_ids |= curr_cell_ids
    trj.reset_index(drop=True, inplace=True)  # Reset indexes, just in case
