#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import sys
import re

def load_sum_rev(sum_file):
    '''
    out.sum file after hypoinverse run
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
        if abs(tmplt_stlo - tar_stlo)<0.04 and abs(tmplt_stla - tar_stla) <0.04:
            s+= f"{line}\n"
    pipe = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE).communicate(s.encode())
    result = re.split(r"\n",pipe[0].decode('utf-8'))
    with open(xc_file,'a') as f:
        for line in result[1:-1]:
            f.write(content[i]+" "+line+"\n")
    f.close()


