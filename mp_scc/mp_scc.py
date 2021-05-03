#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
#     Author: Jinping ZI
#     History: 
#             2021-2-26 Initiate coding
#-----------------------------------------------------------------------------

import os
import subprocess
import multiprocessing as mp
import time
import re
from obspy import UTCDateTime
from utils.basic_utils import str2time,time2str

def load_sum_rev(sum_file):
    '''
    Load in information from out.sum file after hypoinverse run
    format: {eve_folder:eve_id,evlo,evla}
    '''
    sum_list = {}
    with open(sum_file,'r') as f:
        for line in f:
            eve_id=int(line[136:146])
            _e_time = line[0:16]
            e_time = str2time(_e_time)
            e_time = e_time - 8*60*60
            eve_folder = time2str(e_time)[:16]
            evla = int(line[16:18])+0.01*int(line[19:23])/60
            evlo = int(line[23:26])+0.01*int(line[27:32])/60
            sum_list[eve_folder] = [eve_id,evlo,evla]
    f.close()
    return sum_list


def scc_c(sta_pha,point1,point2):
    '''
    SCC written by C is used here because C runs much faster than python. 
    Even though python SCC script has been developed.
    '''
    sum_rev_dict = load_sum_rev("out.sum") # Read in event information
    print("Templates range: %d %d " %(point1,point2))
    xc_file =f"{sta_pha}/{sta_pha}.{format(point1,'05d')}.xc"
    f = open(xc_file,'w')                  # Initiate result file
    f.close()

    # Here we re-read arr files rather than deliver parameter from mp_scc()
    content= []
    with open(f"{sta_pha}.arr","r")as f:
        for line in f:
            line = line.rstrip()
            content.append(line)
    f.close()
    if sta_pha[-1]=="S":
        cmd ="scc -C0.6 -M3 -W-1/3/2 \n"
        #-C0.6:    threshold
        #-M3:      3 components
        #-W-1/3/2: waveform range [-1,3] from arrival time. 
        #          expand sliding range 2s 
        # Run ./SCC in command line to read usage guideline
    if sta_pha[-1]=="P":
        cmd ="scc -C0.6 -M3 -W-0.5/1/0.75 \n"
        # P segment should be short than S waveform
    for i in range(point1,point2):
        s = content[i]+"\n"
        e_path = re.split(" +",line)[0] # Format is "Path Arri_time 1"
        tmplt_folder = re.split("\/",e_path)[-2]
        tmplt_stlo = sum_rev_dict[tmplt_folder][1]
        tmplt_stla = sum_rev_dict[tmplt_folder][2]
        for line in content[i+1:]:
            e_path = re.split(" +",line)[0] # Format is "Path Arri_time 1"
            tar_folder = re.split("\/",e_path)[-2]
            tar_stlo = sum_rev_dict[tar_folder][1]
            tar_stla = sum_rev_dict[tar_folder][2]
            # Accept pairs with 0.04 degree range
            if abs(tmplt_stlo-tar_stlo)<0.04 and abs(tmplt_stla-tar_stla)<0.04:
                s+= f"{line}\n"
        pipe = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE).communicate(s.encode())
        result = re.split(r"\n",pipe[0].decode('utf-8'))
        with open(xc_file,'a') as f:
            for line in result[1:-1]:
                f.write(content[i]+" "+line+"\n")
        f.close()

def mp_scc(sta_pha,cores):
    '''
    The calculation proess is:
        -- The first waveform cross-correlate with left n-1 waveforms
        -- The i-th waeform cross-correlate with left n-j waveforms
    The multiprocessing solution is dividing the whole process into sub-process:
        -- The 1st-200th templates round as one group
        -- The 200th-450th templates round as the second group
        -- The gap increases to maintain general equivalent computation time 
    '''
    s_point = 0          # Start from the first event
    content = []
    with open(f"{sta_pha}.arr","r")as f:
        for line in f:
            line = line.rstrip()
            content.append(line)
    f.close()
    total_amount = len(content)
    gap = 200                      # Will increase in further group
    index_list = [s_point]         # Initiate the list
    loop_j = s_point+gap           # The loop parameter
    while loop_j < total_amount:
        index_list.append(loop_j)
        gap = gap + 50
        loop_j += gap
    index_list.append(total_amount)# Append the end value
    tasks=[]

    for i in range(len(index_list)-1):
        point1 = index_list[i]
        point2 = index_list[i+1]
        tasks.append((sta_pha,point1,point2))
    
    # Multiprocessing
    pool = mp.Pool(processes=cores)
    rs = pool.starmap_async(scc_c,tasks,chunksize=1)
    pool.close()
    while(True):
        remaining = rs._number_left
        print(f"Finished {len(tasks)-remaining}/{len(tasks)}",end = '\r')
        if(rs.ready()):
            break
        time.sleep(0.5)

if __name__ == "__main__":
    """
    Description:
    Scripts for multiprocessing SCC.
    Usage: python mp_scc.py
    The program recognize "*.arr" as arrival files
    """
    af_folder = "./"                         # Arrive files folder
    sta_pha_list = []                        # List for sta_pha
    for file in os.listdir(af_folder):
        if file[-3:]!="arr":                 # Pass non-arrival files
            continue
        sta_pha = re.split("\.",file)[0]     # E.g. sta_pha = GS010_P
        sta_pha_list.append(sta_pha)
    for sta_pha in sta_pha_list:
        os.makedirs(sta_pha)                 # Error happens when exists
        cores = int(mp.cpu_count()/2)        # Expand if you want
        mp_scc(sta_pha,cores)
        finish_time = UTCDateTime.now()
        with open("mp_scc.log",'a') as f:
            f.write(str(finish_time)+" "+sta_pha+"\n")
        f.close()

            

