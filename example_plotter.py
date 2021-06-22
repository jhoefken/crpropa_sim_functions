import matplotlib
#matplotlib.use('Agg') # Uncomment this if you are using Windows
from crpropa import *
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sim_functions as sf
import os

# Distances source
mydata = "distances_v1.dat"
distances = np.genfromtxt(mydata)
distances = distances[:600]  # just considering distances up to z = 0.6
distances = distances.T

# Main parameters
gamma = 1.77
rcut = 20.49
parts = 10

fh = 1.0
fhe = fn = fsi = ffe = 0.

# File names
title = 'output'
plotfile = 'my_plot'
plottitle = 'My simulation and Pierre Auger data'

params = np.array([gamma,rcut,fh,fhe,fn,fsi,ffe])
sf.simulate_dd_1D_parts(params,distances,title=title,parts=parts)
sf.plot_errors_parts(title=title,plotfile=plotfile,plottitle=plottitle,parts=parts)
