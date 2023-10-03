### TSS GC accessibility plots
### Input needs a relative position (either binned or relative to the TSS), mean condition 1, sd condition 1, mean condition 2, sd condition, etc.
### Output is a pretty line graph
### JRBN

# Load libraries
import argparse
import numpy
import matplotlib
import pandas
import math
import os
import csv
import scipy
from scipy import stats, integrate
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon
from matplotlib.table import Table

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generate a TSS GC accessibility plots FANS and 22C 100U")
parser.add_argument('--input', dest="inFile", required=True, help="Input counts CSV")
parser.add_argument('--condition', dest="plotCondit", required=True, help="Condition to plot")
parser.add_argument('--output', dest="outPNG", required=True, help="Output file name base")
args = parser.parse_args()
# Import data


dataIN=pandas.read_csv(args.inFile, header='infer', delimiter=',')

dataINsub=dataIN[['pos', 'mean_'+str(args.plotCondit), 'sd_'+str(args.plotCondit)]].to_numpy()


# Set up data arrays for 22C all
pos_array=[]
mean_array=[]
sd_array=[]
mean_plus_sd_array=[]
mean_minus_sd_array=[]


for cnt in range(0,len(dataINsub)):
    if numpy.isnan(dataINsub[cnt][1]) == False :
        pos_array.append(dataINsub[cnt][0])
        mean_array.append(dataINsub[cnt][1])
        sd_array.append(dataINsub[cnt][2])
        mean_plus_sd_array.append(dataINsub[cnt][1]+dataINsub[cnt][2])
        mean_minus_sd_array.append(dataINsub[cnt][1]-dataINsub[cnt][2])

fig,ax = plt.subplots(figsize=(6,4), dpi=300)
pltLine = ax.plot(pos_array, mean_array, 'b-')
plt.fill_between(pos_array, mean_minus_sd_array, mean_plus_sd_array, alpha=0.25, edgecolor='b', facecolor='b',
                 linewidth=2, linestyle='dashdot', antialiased=True)
ax.grid(color=(0.75,0.75,0.75,0.25), linestyle='-', linewidth=0.5)
plt.xlim=(-1000,1000)
plt.ylim=(-0.1,1.1)
ax.set_ylim(-0.1,1.1)
ax.set_xlabel('Position relative to TSS')
ax.set_ylabel("GC accessibility (%)")
ax.set_title("GC accessibility around TSS (+/- 1kb) -- " + str(args.plotCondit))

#plt.show()

fig.savefig(args.outPNG)


