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
import scipy
from scipy import stats, integrate
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon
from matplotlib.table import Table

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Generate a TSS GC accessibility plots FANS and 22C 100U")
    parser.add_argument('--input', dest="inFile", required=True, help="Input counts CSV")
    parser.add_argument('--output', dest="outPNG", required=True, help="Output file name base")
    args = parser.parse_args()
    return args

def main():
    # Import data
    dataIN=numpy.genfromtxt(args.inFile,delimiter=",", names=True,
                            dtype=[('pos','<f8'),('type','<f8'),('freq','<f8'),('mean_22C_100U','<f8'),('sd_22C_100U','<f8'),
                                   ('mean_FANS_0p5U','<f8'),('sd_FANS_0p5U','<f8'),('mean_FANS_1p5U','<f8'),('sd_FANS_1p5U','<f8'),
                                   ('mean_FANS_25U','<f8'),('sd_FANS_25U','<f8'),('mean_FANS_5U','<f8'),('sd_FANS_5U','<f8')])

    #print(dataIN)

    # Set up data arrays for 22C all
    pos_22C100U_all=[]
    mean_22C100U_all=[]
    sd_22C100U_all=[]
    mean_22C100U_all_plus_sd=[]
    mean_22C100U_all_minus_sd=[]
    
    for cnt in range(0,len(dataIN)):
        if dataIN[cnt][3] != "" :
            pos_22C100U_all.append(dataIN[cnt][0])
            mean_22C100U_all.append(dataIN[cnt][3])
            sd_22C100U_all.append(dataIN[cnt][4])
            mean_22C100U_all_plus_sd.append(dataIN[cnt][3]+(dataIN[cnt][4]))
            mean_22C100U_all_minus_sd.append(dataIN[cnt][3]-(dataIN[cnt][4]))
    
    fig,ax = plt.subplots(figsize=(6,4), dpi=300)
    pltLine22C100Uall=ax.plot(pos_22C100U_all, mean_22C100U_all, 'b-')
    plt.fill_between(pos_22C100U_all, mean_22C100U_all_minus_sd,mean_22C100U_all_plus_sd, alpha=0.25, edgecolor='b', facecolor='b',
                     linewidth=2, linestyle='dashdot', antialiased=True)
    ax.grid(color=(0.75,0.75,0.75,0.25), linestyle='-', linewidth=0.5)
    plt.xlim=(-1000,1000)
    plt.ylim=(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('Position relative to TSS')
    ax.set_ylabel("GC accessibility (%)")
    ax.set_title("GC accessibility around TSS (+/- 1kb) -- 22C 100U GC methylase")
    fig.savefig(args.outPNG + "_22C_100U.pdf")


    # Set up data arrays for 22C all
    pos_22C_all=[]
    mean_22C_all=[]
    sd_22C_all=[]
    mean_22C_all_plus_sd=[]
    mean_22C_all_minus_sd=[]
    
    for cnt in range(0,len(dataIN)):
        if dataIN[cnt][3] != "" :
            pos_22C_all.append(dataIN[cnt][0])
            mean_22C_all.append(dataIN[cnt][5])
            sd_22C_all.append(dataIN[cnt][6])
            mean_22C_all_plus_sd.append(dataIN[cnt][5]+(dataIN[cnt][6]))
            mean_22C_all_minus_sd.append(dataIN[cnt][5]-(dataIN[cnt][6]))
    
    fig,ax = plt.subplots(figsize=(6,4), dpi=300)
    pltLine22Call=ax.plot(pos_22C_all, mean_22C_all, 'b-')
    plt.fill_between(pos_22C_all, mean_22C_all_minus_sd,mean_22C_all_plus_sd, alpha=0.25, edgecolor='b', facecolor='b',
                     linewidth=2, linestyle='dashdot', antialiased=True)
    ax.grid(color=(0.75,0.75,0.75,0.25), linestyle='-', linewidth=0.5)
    plt.xlim=(-1000,1000)
    plt.ylim=(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('Position relative to TSS')
    ax.set_ylabel("GC accessibility (%)")
    ax.set_title("GC accessibility around TSS (+/- 1kb) -- FANS 0.5U")
    fig.savefig(args.outPNG + "_FANS_0p5U.pdf")

    # Set up data arrays for 4C all
    pos_4C_all=[]
    mean_4C_all=[]
    sd_4C_all=[]
    mean_4C_all_plus_sd=[]
    mean_4C_all_minus_sd=[]

    
    for cnt in range(0,len(dataIN)):
        if dataIN[cnt][5] != "" :
            pos_4C_all.append(dataIN[cnt][0])
            mean_4C_all.append(dataIN[cnt][7])
            sd_4C_all.append(dataIN[cnt][8])
            mean_4C_all_plus_sd.append(dataIN[cnt][7]+(dataIN[cnt][8]))
            mean_4C_all_minus_sd.append(dataIN[cnt][7]-(dataIN[cnt][8]))
    
    fig,ax = plt.subplots(figsize=(6,4), dpi=300)
    pltLine4Call=ax.plot(pos_4C_all, mean_4C_all, 'b-')
    plt.fill_between(pos_4C_all, mean_4C_all_minus_sd,mean_4C_all_plus_sd, alpha=0.25, edgecolor='b', facecolor='b',
                     linewidth=2, linestyle='dashdot', antialiased=True)
    ax.grid(color=(0.75,0.75,0.75,0.25), linestyle='-', linewidth=0.5)
    plt.xlim=(-1000,1000)
    plt.ylim=(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('Position relative to TSS')
    ax.set_ylabel("GC accessibility (%)")
    ax.set_title("GC accessibility around TSS (+/- 1kb) -- FANS 1.5U")
    fig.savefig(args.outPNG + "_FANS_1p5U.pdf")


    #print(dataIN)
    # Set up data arrays for 22C common
    pos_22C_common=[]
    mean_22C_common=[]
    sd_22C_common=[]
    mean_22C_common_plus_sd=[]
    mean_22C_common_minus_sd=[]

    
    for cnt in range(0,len(dataIN)):
        if dataIN[cnt][7] != "" :
            pos_22C_common.append(dataIN[cnt][0])
            mean_22C_common.append(dataIN[cnt][9])
            sd_22C_common.append(dataIN[cnt][10])
            mean_22C_common_plus_sd.append(dataIN[cnt][9]+(dataIN[cnt][10]))
            mean_22C_common_minus_sd.append(dataIN[cnt][9]-(dataIN[cnt][10]))
    
    fig,ax = plt.subplots(figsize=(6,4), dpi=300)
    pltLine22Ccommon=ax.plot(pos_22C_common, mean_22C_common, 'b-')
    plt.fill_between(pos_22C_common, mean_22C_common_minus_sd,mean_22C_common_plus_sd, alpha=0.25, edgecolor='b', facecolor='b',
                     linewidth=2, linestyle='dashdot', antialiased=True)
    ax.grid(color=(0.75,0.75,0.75,0.25), linestyle='-', linewidth=0.5)
    plt.xlim=(-1000,1000)
    plt.ylim=(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('Position relative to TSS')
    ax.set_ylabel("GC accessibility (%)")
    ax.set_title("GC accessibility around TSS (+/- 1kb) -- FANS 25U")
    fig.savefig(args.outPNG + "_FANS_25U.pdf")

    # Set up data arrays for 4C common
    pos_4C_common=[]
    mean_4C_common=[]
    sd_4C_common=[]
    mean_4C_common_plus_sd=[]
    mean_4C_common_minus_sd=[]

    
    for cnt in range(0,len(dataIN)):
        if dataIN[cnt][9] != "" :
            pos_4C_common.append(dataIN[cnt][0])
            mean_4C_common.append(dataIN[cnt][11])
            sd_4C_common.append(dataIN[cnt][12])
            mean_4C_common_plus_sd.append(dataIN[cnt][11]+(dataIN[cnt][12]))
            mean_4C_common_minus_sd.append(dataIN[cnt][11]-(dataIN[cnt][12]))
    
    fig,ax = plt.subplots(figsize=(6,4), dpi=300)
    pltLine4Ccommon=ax.plot(pos_4C_common, mean_4C_common, 'b-')

    plt.fill_between(pos_4C_common, mean_4C_common_minus_sd,mean_4C_common_plus_sd, alpha=0.25, edgecolor='b', facecolor='b',
                     linewidth=2, linestyle='dashdot', antialiased=True)
    ax.grid(color=(0.75,0.75,0.75,0.25), linestyle='-', linewidth=0.5)
    plt.xlim=(-1000,1000)
    plt.ylim=(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.set_xlabel('Position relative to TSS')
    ax.set_ylabel("GC accessibility (%)")
    ax.set_title("GC accessibility around TSS (+/- 1kb) -- FANS 5U")
    fig.savefig(args.outPNG + "_FANS_5U.pdf")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()

    # Calling main script
    main()


