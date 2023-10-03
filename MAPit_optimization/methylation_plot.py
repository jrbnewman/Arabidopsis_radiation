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
parser = argparse.ArgumentParser(description="Line plot for TSS methylation")
parser.add_argument('--input', dest="inFile", required=True, help="Input counts CSV")
parser.add_argument('--output', dest="outPDF", required=True, help="Output file name base")
parser.add_argument('--y-min', dest="yMIN", required=True, help="Y-axis minimum")
parser.add_argument('--y-max', dest="yMAX", required=True, help="Y-axis maximum")
parser.add_argument('--title', dest="plotTitle", required=True, help="Plot title")
parser.add_argument('--y-label', dest="yLabel", required=True, help="Y-axis label")
parser.add_argument('--x-label', dest="xLabel", required=True, help="X-axis label")

args = parser.parse_args()
# Import data


dataIN=pandas.read_csv(args.inFile, header='infer', delimiter=',')

#palette
#colors = ["#0033ff","#ff6a00","#000000","#ff00f7","#ff0000","#000000"]
#colors = ["#0033ff","#ff6a00","#000000","#ff00f7","#ff0000","#000000"]
colors = ["#f07705","#3288f0","#f005a2","#000000","#9b00e3","#000000"]



fig,ax = plt.subplots(figsize=(6,6), dpi=300)
num=0
for column in dataIN.drop('pos',axis=1):
    plt.plot(dataIN['pos'], dataIN[column], marker="", color=colors[num], linewidth=2, alpha=1, label=column)
    num+=1

ax.legend( ncol=2, bbox_to_anchor=(0.5, -0.35), fancybox=True, loc='lower center')


#ax.legend(ncol=len("category"), bbox_to_anchor=(0.5, 1.05),fancybox=True,loc='lower center')

#plt.xlim=(-1000,1000)
#plt.ylim=(args.yMIN,args.yMAX)
ax.set_xlim(-1000,1000)
#ax.set_ylim(args.yMIN,args.yMAX)
#ticks=numpy.arange(fl=args.yMIN, args.yMAX, 0.1)
#plt.yticks()
#plt.ylim(float(args.yMIN),float(args.yMAX))
ax.set_ylim(float(args.yMIN),float(args.yMAX))
ax.set_xlabel(args.xLabel)
ax.set_ylabel(args.yLabel)
ax.set_title(args.plotTitle)
plt.tight_layout()

#plt.show()
fig.savefig(args.outPDF)
