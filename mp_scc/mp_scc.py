#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
#     Author: Jinping ZI
#     History: 
#             2021-02-26 Initiate coding
#             2021-06-14 Update parser
#-----------------------------------------------------------------------------

import os
import subprocess
import multiprocessing as mp
import argparse
import time
import re
import warnings
import shutil
from obspy import UTCDateTime

def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--af_folder",
                        default="./",
                        help="The arrival files folder")
    parser.add_argument("--P_scc",
                        default="scc -C0.6 -M3 -W-0.5/1/0.75",
                        help="Command for P arrival cross-correlation")
    parser.add_argument("--S_scc",
                        default="scc -C0.6 -M3 -W-1/3/2",
                        help="Command for S arrival cross-correlation")
    parser.add_argument("--cpu_cores",
                        default=0,
                        type=int,
                        help="0 indicates using all cores")
    args = parser.parse_args()
    return args

def load_sum_rev(sum_file):
    '''
    Load in information from out.sum file after hypoinverse run
    format: {eve_folder:eve_id,evlo,evla}
    '''
    sum_dict = {}
    with open(sum_file,'r') as f:
        for line in f:
            eve_id=int(line[136:146])
            eve_folder = line[:16]
            evla = int(line[16:18])+0.01*int(line[19:23])/60
            evlo = int(line[23:26])+0.01*int(line[27:32])/60
            sum_dict[eve_folder] = [eve_id,evlo,evla]
    f.close()
    return sum_dict


def scc_c(args,sta_pha,point1,point2):
    '''
    SCC written by C is used here because C runs much faster than python. 
    Even though python SCC script has been developed.
    '''
    sum_rev_dict = load_sum_rev("out.sum") # Read in event information
    print("Templates range: %d %d " %(point1,point2))
    xc_file =f"{sta_pha}/{sta_pha}.{format(point1,'05d')}.xc"
    f = open(xc_file,'w')                  # Initiate result file
    f.close()

    # Here we re-read arr files rather than deliver parameters from mp_scc()
    content= []
    with open(f"{sta_pha}.arr","r")as f:
        for line in f:
            line = line.rstrip()
            content.append(line)
    f.close()
    if sta_pha[-1]=="S":
        cmd =args.S_scc+"\n"
        #-C0.6:    threshold
        #-M3:      3 components
        #-W-1/3/2: waveform range [-1,3] from arrival time. 
        #          expand sliding range 2s 
        # Run ./SCC in command line to read usage guideline
    if sta_pha[-1]=="P":
        cmd = args.P_scc+"\n"
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

def mp_scc(args,sta_pha):
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
        tasks.append((args,sta_pha,point1,point2))
    # Multiprocessing
    cores = args.cpu_cores
    if cores==0:
        cores = mp.cpu_count()
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
    Usage: python mp_scc.py  --P_scc="scc -C0.6 -M3 -W-0.5/1/0.75"
    The program recognize "*.arr" as arrival files
    """
    args = read_args()
    sta_pha_list = []                        # List for sta_pha
    for file in os.listdir(args.af_folder):
        if file[-3:]!="arr":                 # Pass non-arrival files
            continue
        sta_pha = re.split("\.",file)[0]     # E.g. sta_pha = GS010_P
        sta_pha_list.append(sta_pha)
    for sta_pha in sta_pha_list:
        if os.path.exists(sta_pha):
            status = 0
            while status != 1:
                in_parameter = input("The folder has existed, remove it? Y/N\n")
                if in_parameter=="Y":
                    shutil.rmtree(sta_pha)
                    status = 1
                elif in_parameter=="N":
                    status = 1
                    os._exit()
                else:
                    print("Please select Y/N...")
        os.makedirs(sta_pha)                 # Error happens when exists
        mp_scc(args,sta_pha)
        finish_time = UTCDateTime.now()
        with open("mp_scc.log",'a') as f:
            f.write(str(finish_time)+" "+sta_pha+"\n")
        f.close()

            

