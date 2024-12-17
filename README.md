# E. coli survival in response to ciprofloxacin antibiotic stress correlates with increased nucleoid length and effective misfolded protein management. 

**George Butler<sup>1,4</sup>, Julia Bos<sup>2,3,4</sup>, Robert H. Austin<sup>3</sup>, Sarah R. Amend<sup>1</sup>, Kenneth J. Pienta<sup>1</sup>**

<sup><sup>1</sup>Cancer Ecology Center, The Brady Urological Institute,  Johns Hopkins School of Medicine, Baltimore, USA</sup>

<sup><sup>2</sup>Unite Plasticite du Genome Bacterien, Institut Pasteur, UMR3525, CNRS, Paris 75015, France</sup>

<sup><sup>3</sup>Department of Physics, Princetone University, Princetone, USA</sup>

<sup><sup>4</sup>Both authors contributed equally</sup>

![quantifying_dispersal_example](/Readme_example_image/segmentation_example.jpg)


## Overview

Code and processed data for the manuscript "Filamentous E. coli revelas interplay of cell responses that increases survival against antibiotic stress", for more information please read our [paper](https://royalsocietypublishing.org/doi/full/10.1098/rsos.230338)

The seperated stages of the segementation pipeline are provided in directories 1_ through 6_ and an example set of raw images are provided in the "demo" directory "example_images.zip" along with an example video. Each stage of the pipeline should be used sequentially and a brief overview of each stage can be found below. The example images can be used from stage 2 onwards in the pipeline at the output from each stage should be used in the next stage e.g the example images should be used as an input in stage 2 and then the output should be used as an input for stage 3 etc. Finally, the directory "statatical analysis" contains the necessary data and R code to reproduce the results in presented in the manuscript. 

Note: The pipeline is designed to run on a UNIX machine and take advantage of multiple threads where possible. All of the stages *should* run on a Windows machine as well but please be aware that the threaded backend will not work due to the Python Global Interpreter Lock, further information can be found [here](https://joblib.readthedocs.io/en/latest/parallel.html)

## Installation
Python 3.6 is required to run the scripts in directories 1_ through 6_. We strongly recommend that you use [Anaconda](https://www.anaconda.com/) to create a virtual environment to prevent conflicts with versions of Python that may be installed on your local system. A virtual environment can be created by running the following from you terminal: 

conda create --name bacteria_filamentation python=3.6

The environment can be activated with:

conda activate bacteria_filamentation 

Finally, the necessary libraries can be install be navigating to the location of the requirements.txt file and running:

pip install -r requirements.txt 
	

## Running each section of the pipeline

1. 1_tiff_seperation

	If the time-lapse data has been collected on a Nikon microscope and saved as an .nd2 file then this script will split the file into seperate images for each channel and XY position. To run the script point to the location of the directory in which the .nd2 file has been saved, referred to from here on as the master directory. Once finished, the seperated images will be stored in a directory called "images" that will be alongside the original .nd2 file. An example of the output can be seen in example_images.zip.
	
2. 2_segemention

	The tiff_to_mask.py script will segment the images that have been split in stage 1. To run the script point to the location of the master directory that contains the images subdirectory and the original .nd2 file. This script can be tested with the example images provided via:
	
	python tiff_to_mask.py ~/bacteria_filamentation/example_images/
	
	The segmented output for each channel / XY position will be saved in a directory called "raw_mask_avg" alongside the original raw images

3. 3_tracking

	The cells in each frame in the phase contrast channel can be tracked by navigating to "ecoli_tracker" directory. The user interface for the tracker can be launched by running:
	
	python cell_main.py
	
	Once the interface has been launched, cells in the phase contrast channel can be tracked by clicking the "Open Folder" button at the top right of the panel. A pop out window will then launch enabling you to navigate to the location of the phase constrast time-lapse images which can be opened by clicking on the "raw" directory and selecting "open" at the bottom right of the pop out window. If the "raw" and "raw_mask_avg" subdirectories are stored with in the same parent directory then the raw phase constrast images will appear on the left of the panel and the segmented images will appear on the right. 
	
	To track the cells click on the "Cell Tracking" button at the upper right of the panel. Once finished, the segmeneted greyscale image will be replaced by a colorised mask where a unique color is chosen for each cell. The individual cells can then be filtered by deselected a given ID on the right hand panel. 
	
	Finally, the output can be saved by clicking the "Save selection" button at the upper right of the panel. To run the further stages in the pipeline the results needed to be saved with the "Show IDs" button toggled on and off.

4. 4_manual_corrections (optional)

	The tracking_postprocessor.py script is optional and designed to account for small unavoidable errors that emerge during the tracking process. The tracking_postprocessing.py script uses a manual correction file as an input called "manual_notes.csv". An **example** manual correation file can be found in the directory. The manual corrections can be determined by opening by the raw and tracked images alongside one another in Fiji as combined image stackes. The key in the manual_notes.csv corresponds to the following situations:
	
	1 = two cells have divided e.g. cell 1 divided in frame 2 into cells 1 and 2 
	
	2 = a cell has been mistracked a different cell e.g. cell 6 should be tracked as cell 3 
	
	3 = a cell needs to be remove from frame X e.g. cell 4 should be removed from frame 5
	
	4 = a cell needs to be removed from frame X and all future frames e.g. cell 6 should be remove from frmae 7 and all subsequent frames 
	
	5 = a cell needs to be removed from all frames e.g. cell 7 should be removed frame all frames (this is the same as deselecting a cell in stage 3). 
	
	Finally, the tracking_postprocessor.py script can be run by pointing the XY directory that contains the raw, raw_mask, labeled, and unlabeled directories. 
										
5. 5_contour_postprocessing

	The contour_postprocessing.py script can be used to extract the shape (or contour) of each tracked cell in every frame. The contour_postprocessing.py script can be run by editing line 186 to record the XY positions that should be used. The script can then be run by pointing to the master directory. 
	
	If you also want to track the indiviudal protein foci then navigate to the "protein_tracking" subdirectory and descend into the "fl_ecoli_tracker" directory and follow the step in stage 3. Note, this is the same tracking process as stage 3 but with different assumptions to account of the smaller protein foci. Once tracked, the protein foci can also be extracted from the images by running the protein_tracker.py script. The protein_tracker.py script functions in the same format was the contour_postprocessing.py script. 
	
6. 6_feature_quantification

	The cell_features.py script can be used to measure the cell area, perimeter and medial_axis. The script can be run by changing lines 56 and 59 to the number of XY positions and the number location of the protein_tracker.py output respectively. 

## Statistical analysis

The "statistical_analysis" directory contains the postprocessed tracking information for each cell as well as the corresponding R script to recreate the Figures and analysis present in our manuscript. 

## Licenses
The copyright of PyQT belong to Riverbank Computing ltd.

Pandas is released under [BSD 3-Clause License](http://pandas.pydata.org/pandas-docs/stable/overview.html?highlight=bsd). Copyright owned by AQR Capital Management, LLC, Lambda Foundry, Inc. and PyData Development Team. 

Trackpy is released under [BSD 3-Clause License](https://github.com/soft-matter/trackpy/blob/master/LICENSE). Copyright owned by trackpy contributors.

NumPy and SciPy are released under BSD-new License

Scikit-image is released under [modified BSD license](https://github.com/scikit-image/scikit-image)

PIMS is released under [modified BSD license](https://github.com/soft-matter/pims/blob/master/license.txt)

Matplotlib is released under [Python Software Foundation (PDF) license](https://matplotlib.org/)

Seaborn is released under [BSD 3-clause license](https://github.com/mwaskom/seaborn/blob/master/LICENSE)

PyQtGraph is released under MIT license
