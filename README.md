# Image trajectories' clustering

This repository contains a group of scripts that allows you to generate and cluster the images of sperms's trajectories. Please follow in order the next steps.

## Generate a *result file* per sample

Firstly, go to the ImageJ folder and install ```BGM_7.java``` as an ImageJ's plug-in. Then, analyze each sperm's samples with this plug-in.
- Input: A sequence of images of a sperm sample.
- Output: A text file (*results file*) that contains the motility parameters and the coordinate vector of the tracked spermatozoa per sample. 

## Group the *result files* into *stacked files*

Run ```log_motility_data-6.py``` for each "results file"  in order to group its data into 6 *stacked files*. Please refer to the folder ```demo_files``` for the standard folder of the input (```results_file.txt```) and output (```Motility_Results.txt```, ```Motility_Results_Avg.txt```, ```Motility_Results_Median.txt```, ```Motility_Results_Sigma.txt``` and ```Motility_Results_Speed.txt```) files.
- Input: All the *results files*. Example for a single results file processing in Ubuntu's terminal ```$ python3 log_motility_data-6.py results_file.txt id1 id2 id3 id4```.
- Output: A *trajectory file* and a *motility parameter file* that contain, respectively, the coordinate vector and the motility parameter of each tracked cell. Additionally, you will obtain 4 text files that contain the average, median, standard deviation and speed of the motility parameters for each sample.

## Generate the *image dataset*

Run ```image_dataset_gen.ipynb``` to generate the *image dataset* that contains the *trajectories' images* extracted from the coordinate vectors stored in the *trajectory file*.
- Input: The *trajectory file*.
- Output: The *image dataset* and the *master dataframe* (optional, but recommended).

## Compute the clusters and update the *master dataframe*

Run ```traj_agglo.ipynb``` in order to group the *trajectories' images* using the hierarchical agglomerative algorithm.
- Input: The *master dataframe* and the number of clusters.
- Output: The *clustered dataframe*, which is the updated version of the *master dateframe* where the clusters are included.

## TODO

- Add links to provide the demo files
- Visualization functions
- Statistical functions for clustered_df
