import pims
import numpy as np
from os.path import join, split, normpath, exists
from os import makedirs
from datetime import datetime
from cell_drawing import create_colorized_masks
from imageio import imwrite
import logging


def read_img_sequence(path, file_extension):
    pims_sequence = pims.ImageSequence(join(path, '*.{}'.format(file_extension)), process_func=None)
    return np.stack([frame.copy() for frame in pims_sequence], axis=2)


def save_results(path_in, trj, col_tuple, col_weights, id_masks, cell_ids, background_id,
                 color_list, cell_color_idx, cell_visibility,
                 show_ids, show_contours,
                 mask_extension):
    # Create output folder based on input path and experiment date/time
    save_time = datetime.now()  # Experiment time
    root_path, folder_in = split(normpath(path_in))
    if show_ids == True:
        path_out = join(root_path,'labeled')
    if show_ids == False:
        path_out = join(root_path, 'unlabeled')

    if not exists(path_out):
        makedirs(path_out)

    save_results_to_csv(path_out, trj, col_tuple, cell_visibility)
    colorized_masks = create_colorized_masks(id_masks, trj, cell_ids, background_id,
                                             color_list, cell_color_idx,
                                             cell_visibility, show_ids)
    print(show_ids)
    logging.info('================= Saving colorized masks frame by frame =================')
    save_sequence_frame_by_frame(colorized_masks, path_out, 'Masks_per_frame', mask_extension, 'masks')
    logging.info('Saving finished.')


def save_results_to_csv(path_out, trj, col_tuple, cell_visibility):
    cols_to_save = ['particle'] + ['frame'] + ['y'] + ['x']
    order_list = ['particle', 'frame']
    trj[trj['particle'].isin([cell_id
                                            for cell_id, show_cell in cell_visibility.items()
                                            if show_cell])].sort_values(order_list).to_csv(
        join(path_out, 'tracks.csv'),
        columns=cols_to_save,
        float_format='%f',  # Use '%.03f' for 3 digits after the comma
        index=False)


def save_sequence_frame_by_frame(sequence, path_out, sequence_folder, file_extension, file_prefix):
    path_out_sequence = join(path_out, sequence_folder)
    if not exists(path_out_sequence):
        makedirs(path_out_sequence)
    n_frames = len(sequence)
    max_n_digits = 1 + int(np.floor(np.log10(n_frames - 1)))
    for i_frame, frame in enumerate(sequence):
        logging.info('Frame {}...'.format(i_frame))
        frame_name = '{}_{}.{}'.format(file_prefix, str(i_frame).zfill(max_n_digits), file_extension)
        imwrite(join(path_out_sequence, frame_name), frame)
