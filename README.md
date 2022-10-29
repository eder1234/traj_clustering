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

## 
