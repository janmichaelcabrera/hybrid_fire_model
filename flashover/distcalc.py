import yaml
import os
import subprocess
import pandas as pd
import numpy as np
import pdb




def main(configlocation,configkey):
    configs = yaml.load(open(configlocation))
    configuse = configs[configkey]
    #iterating over all layouts
    for filename in os.listdir(configuse['layoutloc']):
        if filename.endswith('.csv'):
            rawlayout= pd.read_csv(configuse['layoutloc'] \
            + '/' + filename)
            num_items = len(rawlayout.item)
            distmatrix = np.zeros((num_items,num_items))
            #iterating through all items as sources
            for (i,source) in enumerate(rawlayout.item):
                center = np.array([np.average([rawlayout.iloc[i]['x1']\
                ,rawlayout.iloc[i]['x2']]),np.average([rawlayout.iloc[i]['y1']\
                ,rawlayout.iloc[i]['y2']])])
                #iterating through all j items to determine tendency for 
                #item i to ignite item j
                for (j,item) in enumerate(rawlayout.item):
                    #Finding the x coord closest to source
                    if rawlayout.iloc[j]['x1'] <= center[0] <= rawlayout.iloc[j]['x2']:
                        x_edge = center[0]
                    elif rawlayout.iloc[j]['x1'] > center[0]:
                        x_edge = rawlayout.iloc[j]['x1']
                    elif rawlayout.iloc[j]['x2'] < center[0]:
                        x_edge = rawlayout.iloc[j]['x2']

                    #Finding y coord closest to source
                    if rawlayout.iloc[j]['y1'] <= center[1] <= rawlayout.iloc[j]['y2']:
                        y_edge = center[1]
                    elif rawlayout.iloc[j]['y1'] > center[1]:
                        y_edge = rawlayout.iloc[j]['y1']
                    elif rawlayout.iloc[j]['y2'] < center[1]:
                        y_edge = rawlayout.iloc[j]['y2']
                    distmatrix[i,j] =np.linalg.norm(center-[x_edge,y_edge])
            savename = configuse['layoutloc'] \
            + '/dist' + '/' + filename
            np.savetxt(savename,distmatrix,delimiter=',')


if __name__ == "__main__":
    filename = os.path.splitext(os.path.basename(__file__))[0]
    main('./modelconfig.yaml',filename)