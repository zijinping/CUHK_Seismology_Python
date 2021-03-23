#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
#     Author: Jinping ZI
#     History: 
#             2021-2-26 Initiate coding
#     Usage:
#             e.g. "python mp_scc.py GS010_S"  
#-----------------------------------------------------------------------------

import os
import subprocess
import multiprocessing as mp
import sys
import time

def load_sum_rev(sum_file):
    '''
    Load in information from out.sum file after hypoinverse run
    format: {eve_folder:eve_id,evlo,evla}
    '''
    sum_list = {}
    with open(sum_file,'r') as f:
        for line in f:
            eve_id=int(line[136:146])
            eve_folder = line[0:16]
            evla = int(line[16:18])+0.01*int(line[19:23])/60
            evlo = int(line[23:26])+0.01*int(line[27:32])/60
            sum_list[eve_folder] = [eve_id,evlo,evla]
    f.close()
    return sum_list


def scc_c(sta_pha,point1,point2):
    os.system(f"python sta_scc.py {sta_pha} {point1} {point2}")
    sum_rev_dict = load_sum_rev("out.sum")
    sta_pha = sys.argv[1]
    point1 = int(sys.argv[2])
    point2 = int(sys.argv[3])
    print(point1,point2)
    xc_file =f"{sta_pha}.{format(point1,'05d')}.xc"
    f = open(xc_file,'w')
    f.close()

    content= []
    with open(f"{sta_pha}.arr","r")as f:
        for line in f:
            line = line.rstrip()
            content.append(line)
    if sta_pha[-1]=="S":
        cmd ="scc -C0.6 -M3 -W-1/3/2 \n"
    if sta_pha[-1]=="P":
        cmd ="scc -C0.6 -M3 -W-0.5/1/0.75 \n"
    for i in range(point1,point2):
        s = content[i]+"\n"
        tmplt_folder = content[i][21:37]
        tmplt_stlo = sum_rev_dict[tmplt_folder][1]
        tmplt_stla = sum_rev_dict[tmplt_folder][2]
        for line in content[i+1:]:
            tar_folder = line[21:37]
            tar_stlo = sum_rev_dict[tar_folder][1]
            tar_stla = sum_rev_dict[tar_folder][2]
            if abs(tmplt_stlo-tar_stlo)<0.04 and abs(tmplt_stla-tar_stla)<0.04:
                s+= f"{line}\n"
            pipe = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE).communicate(s.encode())
            result = re.split(r"\n",pipe[0].decode('utf-8'))
            with open(xc_file,'a') as f:
                for line in result[1:-1]:
                    f.write(content[i]+" "+line+"\n")
            f.close()

def mp_scc(sta_pha,cores)
    start_point = 0
    
    content = []
    with open(f"{sta_pha}.arr","r")as f:
        for line in f:
            line = line.rstrip()
            content.append(line)
    f.close()
    total_amount = len(content)
    gap = 200
    index_list=[start_point,start_point+gap]
    sum=start_point+gap
    while sum < total_amount-1:
        gap = gap + 50
        sum+=gap
        if sum <=total_amount:
            index_list.append(sum)
    index_list.append(total_amount)

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
    Usage: python mp_scc.py {sta}_{phase}. Phase could be eight P or S.
    There should be an arr_files folder with station arrive file named in
        {sta}_{phase}.arr
    The .arr file should contain 4 columns:
        1.

    """
    af_folder = "arr_files"                  #Arrive files folder
    sta_pha_list = []                        #List for sta_pha
    for file in os.listdir(af_folder):
        if file[-3:]!=".arr":                #Arrive files ".arr" suffix
            continue
        sta_pha = re.split("\.",file)[0]     #E.g. GS010_P.arr
        sta_pha_list.append(sta_pha)
    for sta_pha in sta_pha_list:
        os.makedirs(sta_pha)                 #Error happen when exists

    sta_pha = sys.argv[1]
    cores = mp.cpu_count()
    mp_scc(sta_pha,cores)
