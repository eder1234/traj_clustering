#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 13:02:02 2022

@author: rodriguez
"""

from matplotlib import pyplot as plt
import pickle

class Trajectory:
    def __init__(self, name, date, quantity, exposure, track_id, traj, 
                 VCL, VAP, VSL, LIN, STR, WOB, BeatCross, ALH, img, 
                 agglo_cluster, SCAN_cluster, SCAN_values, img_size = [144, 144]):
        self.name = name
        self.date = date
        self.quantity = quantity
        self.exposure = exposure
        self.track_id = track_id
        self.traj = traj
        self.VCL = VCL
        self.VAP = VAP
        self.VSL = VSL
        self.LIN = LIN
        self.STR = STR
        self.WOB = WOB
        self.BeatCross = BeatCross
        self.ALH = ALH
        self.img = img
        self.agglo_cluster = agglo_cluster
        self.SCAN_cluster = SCAN_cluster
        self.SCAN_values = SCAN_values
        self.img_size = img_size
        
    def image_display(self, img_size = [144,144]):
        img_reshaped = self.img.reshape(self.img_size)
        plt.imshow(img_reshaped, interpolation='nearest')
        plt.show()
        plt.close()
        
class Dataset:
    def __init__(self, namefile = 'objects/dataset.obj'):
        self.trajectory = []
        self.namefile = namefile
        
    def add_particle(self, trajectory):
        self.trajectory.append(trajectory)
        
    def save(self, namefile = 'objects/dataset.obj'):
        with open(namefile, 'wb') as output:
            pickle.dump(self, output)
            
    def load(self, namefile = 'objects/dataset.obj'):
        with open(namefile, 'rb') as output:
            self = pickle.load(output)
            return self
