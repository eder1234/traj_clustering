# Image trajectories' clustering

This repository contains a group of scripts that allows you to generate and cluster the images of sperms's trajectories. Please follow in order the next steps.

## Generate a *result file* per sample

Firstly, go to the ImageJ folder and install ```BGM_7.java``` as an ImageJ's plug-in. Then, analyze each sperm's samples with this plug-in.
- Input: A sequence of images of a sperm sample.
- Output: A text file (*results file*) that contains the motility parameters and the coordinate vector of the tracked spermatozoa per sample. 

## Group the *result files* into *stacked files*

Run ```log_motility_data.py``` for each "result file"  in order to group its data into 6 *stacked files*.
- Input: All the *results files*.
- Output: A *trajectory file* and a *motility parameter file* that contain, respectively, the coordinate vector and the motility parameter of each tracked cell. Additionally, you will obtain 4 text files that contain the average, median, standard deviation and speed of the motility parameters for each sample.

## Generate the *image dataset*

Run ```image_dataset_gen.ipynb``` to generate the *image dataset* that contains the *trajectories' images* extracted from the coordinate vectors stored in the *trajectory file*.
- Input: The *trajectory file*.
- Output: The *image dataset* and the *master dataframe* (optional, but recommended).

## Compute the clusters and update the *master dataframe* (Pending...)

Run ```traj_agglo.ipynb``` in order to group the *trajectories' images* using the hierarchical agglomerative algorithm.
- Input: ?
- Output: The *master dataframe* updated (clusters included).
