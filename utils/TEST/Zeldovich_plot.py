# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 00:30:40 2025

@author: Tomas
"""
import numpy as np
import h5py as h5
import h5py as h5
pf = h5.File(r"/home/tomas/Desktop/HEDGEHOG_6_9_26/outputDir/3_particles.h5", 'r')
pos_x=pf["pos_x"][:]
pos_y=pf["pos_y"][:]
pos_z=pf["pos_z"][:]
vel_x=pf["vel_x"][:]
vel_y=pf["vel_y"][:]
vel_z=pf["vel_z"][:]
density=pf["density"][:]
current_a = pf.attrs['current_a']
from matplotlib import pyplot as plt
plt.figure(figsize = [10,10])
plt.plot(pos_x[0:64],vel_x[0:64])
plt.scatter(pos_x[0:64],vel_x[0:64])
#plt.plot(pos_x[0:32],density[:,16,16])
#plt.yscale('log')
pf.close()