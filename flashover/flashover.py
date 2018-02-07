import yaml
import os
import subprocess
import pandas as pd
import numpy as np
import pdb

def firecalc(configlocation, room, layout, firstign):
    #Setting the location of config file, which holds settings
    configs = yaml.load(open(configlocation))
    #In that file, the 'firecalc' section is relevant for this code
    configuse = configs['firecalc']
    #Creating an array of times spaced by timestep that goes to simtime
    timelist = np.arange(0,configuse['simtime'], configuse['timestep'])
    #Importing item info as pandas df
    iteminfo = pd.read_csv(configuse['itemloc'],index_col='item')
    #Importing room specs as pandas df
    roominfo = pd.read_csv(configuse['roomloc'] \
    + '/' + room + '.csv')

    

    #Need room surface area for flashover correlation
    surf_area = 2*roominfo['width'][0]*roominfo['length'][0] \
    + 2*roominfo['width'][0]*roominfo['height'][0] \
    + 2*roominfo['length'][0]*roominfo['height'][0]


    #Now calculate the required HRR that produces flashover
    Q_FO = 378*(1+.021*(surf_area/(roominfo['ventarea'][0]*\
        roominfo['ventheight'][0]**.5)))\
        *roominfo['ventarea'][0]*roominfo['ventheight'][0]**.5

    #If/when flashover occurs, the NaN is replaced with a time
    flashover = float('NaN')


    #construct an HRR dictionary holding curves for each item
    hrrdic = {}

    #Incluedes a draw for multiple curves if applicable
    for itemtype in iteminfo.index:
        hrrdic[itemtype] = np.loadtxt(configuse['itemhrrdataloc']\
        +itemtype+ str(np.random.randint(0,iteminfo.numcurves[itemtype])) \
        + '.csv',delimiter=',')

    

    #Retrieving a list of items in the room
    rawlayout= pd.read_csv(configuse['layoutloc'] \
    + '/' + layout + '.csv')
    itemlist = rawlayout.item

    #Retriving the centroid to edge distance matrix for point source
    distmatrix = np.loadtxt(configuse['layoutloc'] \
    + 'dist/' + layout + '.csv', delimiter=',')

    num_items = len(itemlist)

    #Firelist holds the time at which each item is ignited
    firelist = np.ones(num_items)*configuse['simtime']
    firelist[firstign] = 0

    #Creating a ist of accumulated FTP for each item
    FTP = np.zeros(num_items) #FTP for each item


    for t, time in enumerate(timelist):

        #The list of incident fluxes. Resets each time step
        incident_flux = np.zeros(num_items)
        #The current HRR of the total fire, resets at each time step
        HRR = 0

        for f,fire in enumerate(itemlist):
            #Add all fires' contribution to overall HRR
            #Negative times (unignited) are given zero in interp
            HRR = HRR + np.interp(time - firelist[f]\
            ,hrrdic[fire][:,0],hrrdic[fire][:,1])



            #Check if flashover conditions are met
            if HRR > Q_FO and np.isnan(flashover):
                flashover = time
                #If anything hasn't ignited, it has now
                firelist[firelist==configuse['simtime']] = time

            for i,item in enumerate(itemlist):
                #Only executes if ith item hasn't ignited 
                #AND fth item has ignited
                if item == 'HF':
                        incident_flux[i] = incident_flux[i] \
                        + iteminfo.radfrac[fire]\
                        *np.interp(time - firelist[f]\
                        ,hrrdic[fire][:,0],hrrdic[fire][:,1]) \
                        /(4*np.pi*distmatrix[f][i]**2)
                        print(incident_flux[i], i)
                if firelist[f] != configuse['simtime'] and\
                firelist[i] == configuse['simtime']:

                    #Use point source to find incident flux on ith item
                    incident_flux[i] = incident_flux[i] \
                    + iteminfo.radfrac[fire]\
                    *np.interp(time - firelist[f]\
                    ,hrrdic[fire][:,0],hrrdic[fire][:,1]) \
                    /(4*np.pi*distmatrix[f][i]**2)


                    #Then use FTP model
                    if incident_flux[i] > iteminfo.qcrit[item]:
                        FTP[i] = FTP[i] + (incident_flux[i]\
                        -iteminfo.qcrit[item])\
                        **iteminfo.n[item]*configuse['timestep']
             
                        if FTP[i] > iteminfo.FTP[item]:
                            firelist[i] = time
    return [flashover,firelist]

    
                                           
