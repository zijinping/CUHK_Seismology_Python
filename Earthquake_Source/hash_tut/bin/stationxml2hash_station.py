#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stationxml to HASH SCRN.station
    Append station metadata to the desire output file
Created on Thu Mar  5 22:11:49 2020

@author: jw


Old SCEDC format (examples 1-3, 5; station_subs.f; scsn.stations) 
---------------------------------------------------------------------
columns format  value
---------------------------------------------------------------------
1-4     a4      station name
6-8     a3      station component
42-50   f9.5    lattitude (degrees, signed)
52-61   f10.5   longitude (degrees, signed)
63-67   i5      elevation (meters)
91-92   a2      network code
---------------------------------------------------------------------
e.g.ADO  BHE Adelanto Receiving Station       34.55046 -117.43391   878 2000/08/09 3000/01/01 CI

Convert dataless to stationxml:
    https://github.com/iris-edu/stationxml-seed-converter


# MO dataless is broken, thus the station dts manuely add in the station list
DTS  EH* Macou Station                    22.15950  113.56830   134 2000/01/01            MO
"""
from obspy import read_inventory
import os

# Update the desire network
network = "HK"
path_to_stationxml = '/Users/jeremy/OneDrive - The Chinese University of Hong Kong/cu/academic/sources/program/seismology/tmp/Dataless'

outpath = ''
outfile = 'HK'      # without suffix

inv = read_inventory(os.path.join(path_to_stationxml, network + '.xml'))

with open(os.path.join(outpath, outfile + ".stations"), "a") as ostations, \
     open(os.path.join(outpath, outfile + ".reverse"), "a") as oreverse:
    for net in inv:
        # Skip non-selected network 
        if net.code == network:
            pass
        else:
            continue

        for sta in net:
            for cha in sta:
                # Skip specific channels
                if 'HH' in cha.code:
                    pass
                elif 'BH' in cha.code:
                    pass
                else:
                    continue # Skip station not in HH and BH

                # Create station format
                _sta = '{:<4}'.format(sta.code)[:4]
                code = '{:<3}'.format(cha.code)
                try:
                    desc = '{:<33}'.format(' ' + sta.site.name)
                except TypeError:       # Skip empty site name
                    desc = '{:<33}'.format(' ')
                lat  = '{:>9}'.format('%.5f' % cha.latitude)
                lon  = '{:>10}'.format('%.5f' % cha.longitude)
                ele  = '{:>5}'.format('%5d' % cha.elevation)
                try :
                    sdate = '{:<23}'.format(' ' + sta.start_date.strftime('%Y/%m/%d')
                                            + ' ' + sta.end_date.strftime('%Y/%m/%d'))
                except AttributeError:  # Skip empty end date in that station
                    sdate = '{:<23}'.format(' ' + sta.start_date.strftime('%Y/%m/%d'))
                _net = '{:<2}'.format(net.code)

                line = _sta + ' ' + code + desc + lat + ' ' + lon + ' ' \
                       + ele + sdate + _net
                ostations.write(line + '\n')
            oreverse.write('{:<4}'.format(sta.code)[0:4] + ' 0        0 \n')
