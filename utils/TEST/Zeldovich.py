# -*- coding: utf-8 -*-
"""
Created on Wed May  8 23:46:20 2026

@author: Tomas Ondro
"""
# The Zeldovich pancake test
# Provides a good simultaneous test of particle dynamics,
# Poisson solver, and cosmological expansion.
# Valinia et al(1997) Gravitational instability in collisionless cosmological pancakes

import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

###############################################################################
outputFileName = ( '/home/tomas/Desktop/HEDGEHOG_6_9_26/ICs/0_particles.h5' )
Lbox = 64.0          # Mpc/h
n = 128               # particles per dimension
z_c = 9.0            # collapse redshift
z = 19.0             # initial redshift

################################# Cosmology ###################################
H0 = 50.0            # km/s/Mpc
h = H0 / 100.0
Om_M = 1.0
Om_L = 0.0
G = 4.300927161e-9

###############################################################################
a = 1.0 / (1.0 + z)
a_zc = 1.0 / (1.0 + z_c)
rho_crit = ( 3.0 * (H0/h)**2 / (8.0 * np.pi * G) )
rho_cdm_mean = Om_M * rho_crit
dx = Lbox / n
q = ( np.arange(n, dtype=float) * dx + 0.5 * dx )
k = 2.0 * np.pi / Lbox
D = a / a_zc
displacement = np.sin(k*q) / k
x = q + D * displacement
H_a = H0 * np.sqrt( Om_M / a**3 + Om_L )
v = ( a * H_a * D * displacement / h )
N_particles = n**3
pos_x = np.empty(N_particles, dtype=np.float64)
pos_y = np.empty(N_particles, dtype=np.float64)
pos_z = np.empty(N_particles, dtype=np.float64)
vel_x = np.empty(N_particles, dtype=np.float64)
vel_y = np.zeros(N_particles, dtype=np.float64)
vel_z = np.zeros(N_particles, dtype=np.float64)
index = 0
for iy in range(n):
    y = iy * dx + 0.5 * dx
    for iz in range(n):
        zpos = iz * dx + 0.5 * dx
        for ix in range(n):
            pos_x[index] = x[ix]
            pos_y[index] = y
            pos_z[index] = zpos
            vel_x[index] = v[ix]
            index += 1
particle_mass_value = ( rho_cdm_mean * Lbox**3 / N_particles )
particle_mass = np.full( N_particles, particle_mass_value, dtype=np.float64 )
pos_x *= 1000.0
pos_y *= 1000.0
pos_z *= 1000.0
Lbox_kpc = Lbox * 1000.0
pos_x = np.mod(pos_x, Lbox_kpc)
pos_y = np.mod(pos_y, Lbox_kpc)
pos_z = np.mod(pos_z, Lbox_kpc)
print("Number of particles =", N_particles)
print("Particle mass       =", particle_mass_value, "Msun/h")
print("Mean density        =", rho_cdm_mean,
      "(Msun/h)/(Mpc/h)^3")

print("Position range X =", pos_x.min(), pos_x.max(), "kpc/h")
print("Velocity range X =", vel_x.min(), vel_x.max(), "km/s")

with h5.File(outputFileName, 'w') as outFile:
    outFile.create_dataset( 'current_z', data=z )
    outFile.create_dataset( 'current_a', data=a )
    outFile.create_dataset( 'particle_mass', data=particle_mass )
    outFile.create_dataset( 'pos_x', data=pos_x )
    outFile.create_dataset( 'pos_y', data=pos_y )
    outFile.create_dataset( 'pos_z', data=pos_z )
    outFile.create_dataset( 'vel_x', data=vel_x )
    outFile.create_dataset( 'vel_y', data=vel_y )
    outFile.create_dataset( 'vel_z', data=vel_z )