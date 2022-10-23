#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 12:41:00 2022

@author: rodriguez
"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from PIL import Image
from tqdm import tqdm
from numpy import asarray

def master_df(traj_path, mot_path, gen_img=False):
    
    t_df = traj_df(traj_path)
    mr_df = mot_df(mot_path)
    t_df = t_df.sort_values(by=['name', 'date', 'quantity', 'exposure', 'tracked_id']) 
    mr_df = mr_df.sort_values(by=['ID1', 'ID2', 'ID3', 'ID4', 'id'])
    
    traj_name_list = t_df["name"].values.tolist()
    traj_date_list = t_df["date"].values.tolist()
    traj_quantity_list = t_df["quantity"].values.tolist()
    traj_exposure_list = t_df["exposure"].values.tolist()
    traj_tracked_id_list = t_df["tracked_id"].values.tolist()
    traj_traj_list = t_df["traj"].values.tolist()
    mr_VCL_list = mr_df["VCL"].values.tolist()
    mr_VAP_list = mr_df["VAP"].values.tolist()
    mr_VSL_list = mr_df["VSL"].values.tolist()
    mr_LIN_list = mr_df["LIN"].values.tolist()
    mr_STR_list = mr_df["STR"].values.tolist()
    mr_WOB_list = mr_df["WOB"].values.tolist()
    mr_BeatCross_list = mr_df["BeatCross"].values.tolist()
    mr_ALH_list = mr_df["ALH"].values.tolist()
    
    master_dict = {"name":traj_name_list, "date":traj_date_list, "quantity":traj_quantity_list, 
               "exposure":traj_exposure_list, "tracked_id":traj_tracked_id_list, "traj":traj_traj_list,
              "VCL":mr_VCL_list, "VAP":mr_VAP_list, "VSL":mr_VSL_list, "LIN":mr_LIN_list, "STR":mr_STR_list,
              "WOB":mr_WOB_list, "BeatCross":mr_BeatCross_list, "ALH":mr_ALH_list}
    m_df = pd.DataFrame(master_dict)

    if gen_img:            
        master = img_to_folder(m_df)
        master.to_pickle("objects/master.pkl")
    else:
         master = pd.read_pickle("objects/master.pkl")   
        
    return m_df

def traj_df(traj_path):
    df = pd.read_csv(traj_path, names=["name", "date", "quantity", "exposure", "tracked_id", "x", "y"], engine='python') 
    t_df = pd.DataFrame({"name":[], "date":[], "quantity":[], "exposure":[], "tracked_id":[], "traj":[]})

    df1 = df.groupby('name')
    for k1, v1 in df1:
        df2 = df1.get_group(k1)
        df3 = df2.groupby('date')
        for k3,v3 in df3:
            df4 = df3.get_group(k3)
            df5 = df4.groupby('quantity')
            for k5,v5 in df5:
                df6= df5.get_group(k5)
                df7 = df6.groupby('exposure') 
                for k7,v7 in df7:
                    df8 = df7.get_group(k7)
                    df9 = df8.groupby('tracked_id')
                    for k9,v9 in df9:
                        df10 = df9.get_group(k9)
                        x = df10['x'].to_numpy()
                        y = df10['y'].to_numpy()
                        temp_df = pd.DataFrame({"name":[k1], "date":[k3], "quantity":[k5], "exposure":[k7], "tracked_id":[k9], "traj":[[x,y]]})
                        t_df = pd.concat([temp_df, t_df], ignore_index=True)
    return t_df

def mot_df(mot_path):
    mr_df = pd.read_csv(mot_path, engine='python') 
    id_ = [x for x in range(len(mr_df))]
    mr_df["id"] = id_
    return mr_df

def img_to_folder(m_df):
    traj_list = m_df['traj'].values.tolist()
    window_size = 60
    list_arrays = []
    for index, traj in tqdm(enumerate(traj_list)):
        list_x, list_y = zip(traj)
        x = np.asarray(list_x).squeeze()
        y = np.asarray(list_y).squeeze()
        #c = range(len(x)) #arange numpy
        x_c, y_c = center_of_traj(x, y)
        plt.figure(figsize=[2,2])
        plt.plot(x, y, 'k')
        plt.xlim([x_c - window_size/2, x_c + window_size/2])
        plt.ylim([y_c - window_size/2, y_c + window_size/2])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis('off')
        new_id = '{:05d}'.format(index)
        plt.savefig(f'img_traj/{new_id}.png')
        plt.close()
        img = Image.open(f'img_traj/{new_id}.png').convert('L')
        numpydata = asarray(img).flatten()/255
        list_arrays.append(numpydata)
    m_df['img'] = list_arrays
    return m_df

def center_of_traj(x, y):
    max_x = np.amax(x)
    min_x = np.amin(x)
    max_y = np.amax(y)
    min_y = np.amin(y)
    x_c = np.mean([max_x, min_x])
    y_c = np.mean([max_y, min_y])
    return x_c, y_c 

traj_path = "traj (copy).txt"
mot_path = "Motility_Results.txt"

m_df = master_df(traj_path, mot_path)
