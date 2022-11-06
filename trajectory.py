#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 13:02:02 2022

@author: rodriguez
"""

import pickle

class Particle:
    def __init__(self, name, date, quantity, exposure, track_id, traj, 
                 VCL, VAP, VSL, LIN, STR, WOB, BeatCross, ALH, img, 
                 agglo_cluster, SCAN_cluster, SCAN_values):
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
        
class Dataset:
    def __init__(self):
        self.particles = []
        
    def add_particle(self, particle):
        self.particles.append(particle)
        
    def save(self, namefile = 'objects/dataset.obj'):
        with open(namefile, 'wb') as output:
            pickle.dump(self, output)
            
    def load(self, namefile = 'objects/dataset.obj'):
        with open(namefile, 'rb') as output:
            self = pickle.load(output)
            return self