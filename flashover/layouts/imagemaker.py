
import os
import subprocess
import pandas as pd
import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#Specified by user
roomwidth = 6.0 #xdir
roomlength = 5.0 #ydir

#Loops through all csv files in current directory and creates a pictures in the ./pictures directory
for filename in os.listdir('.'):
    if filename.endswith('.csv'):
        layout = pd.read_csv(filename)

        fig1 = plt.figure()

        for num,item in enumerate(layout['item'].values):
            ax1 = fig1.add_subplot(111, aspect='equal')
            ax1.add_patch(
                patches.Rectangle(
                (float(layout['x1'][num]), float(layout['y1'][num])),   # (x,y)
                float(layout['x2'][num])-float(layout['x1'][num]),  # width (x)
                float(layout['y2'][num])-float(layout['y1'][num]),  #length (y)
                fill=False
                )
            )
            centerx = (float(layout['x1'][num])+float(layout['x2'][num]))/2
            centery = (float(layout['y1'][num])+float(layout['y2'][num]))/2
            plt.text(centerx, centery,layout['item'][num])



        plt.xlim(0,roomwidth)
        plt.ylim(0,roomlength)
        plt.xlabel('Width (m)', fontsize=12)
        plt.ylabel('Length (m)', fontsize=12)
        fig1.savefig('./pictures/' + filename[:-4] + '.pdf', dpi=90, bbox_inches='tight')
        plt.close()

