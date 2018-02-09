import yaml
import os
import subprocess
import pandas as pd
import numpy as np
import pdb

def firecalc(configlocation, room, layout, firstign):
    #config location file set 
    configs = yaml.load(open(configlocation))
    #from configs in yaml use 'firecalc' section
    config_use = configs['firecalc']
    #create array of times
    timelist = np.arange(0, config_use['simtime'], config_use['timestep'])
    #import item info as pandas df
    iteminfo = pd.read_csv(config_use['itemloc'], index_col = 'item')
    #import room specs as pandas df
    roominfo = pd.read_csv(config_use['roomloc'] + '/' + room + '.csv')

    #room surface area for flashover correlation
    surf_area = 2*roominfo['width'][0]*roominfo['length'][0] + 2*roominfo['width'][0]*roominfo['height'][0] + 2*roominfo['length'][0]*roominfo['height'][0]

    #find flash over for given compartmetn with surface area = surf_area
    Q_FO = 378*(1 + .021*(surf_area/(roominfo['ventarea'][0] * roominfo['ventheight'][0]**.5))) * roominfo['ventarea'][0]*roominfo['ventheight'][0]**.5

    #if/when flashover occurs, the NaN is replaced with a time
    flashover = float('NaN')

    #retrieve a list of items in the room from the yaml
    raw_layout = pd.read_csv(config_use['layoutloc'] + '/' + layout + '.csv')
    item_list = raw_layout.item

    #finds number of heat flux gauges present
    nHF = 0 
    for n in range(len(item_list)):
        if item_list[n] == 'HF':
            nHF += 1

    #construct an HRR dictionary for room items
    hrr_dic = {}

    #Includes a draw for multiple curves if applicable
    # #the below code as an idea for improving dictionary setup
    # for f, fire in enumerate(item_list):
    #     hrr_dic[fire] = np.loadtxt(config_use['itemhrrdataloc'] + fire + str(np.random.randint(0,iteminfo.numcurves[fire])) + '.csv',delimiter=',')

    for itemtype in iteminfo.index:
        hrr_dic[itemtype] = np.loadtxt(config_use['itemhrrdataloc'] + itemtype + str(np.random.randint(0,iteminfo.numcurves[itemtype])) + '.csv',delimiter=',')

    dist_matrix = np.loadtxt(config_use['layoutloc'] + 'dist/' + layout + '.csv', delimiter = ',')

    num_items = len(item_list)

    #Firelist holds the time at which each item is ignited
    firelist = np.ones(num_items)*config_use['simtime']
    firelist[firstign] = 0

    #Creating a ist of accumulated FTP for each item
    FTP = np.zeros(num_items) #FTP for each item

    heat_flux_array = np.zeros((len(timelist), nHF), dtype = float)

    for t, time in enumerate(timelist):

        #list of incident fluxes, resets at each time step
        incident_flux = np.zeros(num_items)

        #the current HRR of the total fire, resets at each time step
        HRR = 0

        for f, fire in enumerate(item_list):
            # print(fire)
            #add all fire contibutions to overall HRR
            #negative times (unignited) are given zero in interpolation
            HRR = HRR + np.interp(time - firelist[f], hrr_dic[fire][:,0], hrr_dic[fire][:,1])

            if HRR > Q_FO and np.isnan(flashover):
                flashover = time
                #ignite any un-ignited items in compartment
                firelist[firelist==config_use['simtime']] = time

            j=0
            for i, item in enumerate(item_list):
                #calculates the incident heat flux to a particular sensor from all items in room
                if item == 'HF' and fire != 'HF':
                    heat_flux_array[t, j] = iteminfo.radfrac[fire]*np.interp(time - firelist[f],hrr_dic[fire][:,0],hrr_dic[fire][:,1])/(4*np.pi*dist_matrix[f][i]**2)
                    j += 1

                #only executes if ith item has not ignited and fth item has
                if firelist[f] != config_use['simtime'] and firelist[i] == config_use['simtime']:
                    #use point source model to find incident heat flux
                    incident_flux[i] = incident_flux[i] + iteminfo.radfrac[fire]*np.interp(time - firelist[f],hrr_dic[fire][:,0],hrr_dic[fire][:,1])/(4*np.pi*dist_matrix[f][i]**2)

                    #now use Flux Time Product model
                    if incident_flux[i] > iteminfo.qcrit[item]:
                        FTP[i] = FTP[i] + (incident_flux[i] - iteminfo.qcrit[item])**iteminfo.n[item]*config_use['timestep']

                        if FTP[i] > iteminfo.FTP[item]:
                            firelist[i] = time

    return[flashover, firelist, timelist, heat_flux_array]


flashover, firelist, timelist, heat_flux_array = firecalc('modelconfig.yaml', 'burn_structure', 'burn_structure', 0)

