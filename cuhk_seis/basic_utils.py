#/usr/bin/bash
# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
#     Author: Jinping ZI
#
# Revision History
#     2021-01-24 Initiate coding
#-----------------------------------------------------------------------------

import os
import numpy as np
from math import radians,cos,acos,sin,asin,sqrt,ceil,pi
import obspy
from obspy import Stream
import glob
import re
from obspy import UTCDateTime
import pandas as pd

root_path="/NAS1/Sichuan_data/continous_waveform_sac/"

def read_sac_ref_time(tr):
    """
    Read and return reference time of a sac file in obspy.UTCDateTime format.

    Parameter
    --------
    tr: Trace object of obspy
    """

    nzyear = tr.stats.sac.nzyear
    nzjday = tr.stats.sac.nzjday
    nzhour = tr.stats.sac.nzhour
    nzmin = tr.stats.sac.nzmin
    nzsec = tr.stats.sac.nzsec
    nzmsec = tr.stats.sac.nzmsec*0.001
    year,month,day = month_day(nzyear,nzjday)
    sac_ref_time = UTCDateTime(year,month,day,nzhour,nzmin,nzsec)+nzmsec
    return sac_ref_time


def spherical_dist(lon_1,lat_1,lon_2,lat_2):
    """
    Calculate the distance of two postions and return distance in degree
    """

    lon_1,lat_1,lon_2,lat_2 = map(radians,[lon_1,lat_1,lon_2,lat_2])
    a=acos(sin(lat_1)*sin(lat_2)+cos(lat_1)*cos(lat_2)*cos(lon_2-lon_1))
    return a*180/pi


def get_st(net,sta,starttime,endtime,f_folder):
    """
    Read and return waveform between starttime and endtime by specified
    net and station in designated folder. It will merge waveform if include
    more than one file.
    
    The return is a obspy Stream object
    """
    inc_list=[]
    for file in os.listdir(f_folder):
        file_path = os.path.join(f_folder,file)
        try:
            st = obspy.read(file_path,headonly=True)
        except:
            continue
        t1,t2 = st[0].stats.starttime,st[0].stats.endtime
        if t2 < starttime or t1 > endtime \
            or st[0].stats.network != net\
            or st[0].stats.station != sta:
            continue
        else:
            inc_list.append(file_path)
    #Read in data
    st = Stream()
    for path in inc_list:
        st += obspy.read(path)
    if len(st) == 0:
        pass
    else:
        st.trim(starttime,endtime)
    return st

def julday(year,month,day):
    ref_time=UTCDateTime(year,1,1)
    tar_time=UTCDateTime(year,month,day)
    julday=(tar_time-ref_time)/(24*60*60)+1
    return int(julday)

def month_day(year,julday):
    """
    Transfer from julday to month and day.
    Return year,month,day
    """
    #check if year is leap year
    leap=False
    if year%100==0:
        if year%400==0:
            leap=True
    else:
        if year%4==0:
            leap=True
    normal_list=[0,31,59,90,120,151,181,212,243,273,304,334,365]
    leap_list=[0,31,60,91,121,152,182,213,244,274,305,335,366]
    if leap:
        i=0
        while leap_list[i]<julday:
            i=i+1
        month=i
        day=julday-leap_list[i-1]
        return year,month,day
    else:
        i=0
        while normal_list[i]<julday:
            i=i+1
        month=i
        day=julday-normal_list[i-1]
        return year,month,day

def find_nearest(array,value):
    """
    find the nearest value. The return is index and diff
    """
    if type(array) != np.ndarray:
        array=np.array(array)
    idx=np.abs(array-value).argmin()
    diff=array[idx]-value
    return idx,diff

class WY_para():
    '''
    This class reads parameters in "wy.para", which contains parameters for GMT plot.
    '''
    def __init__(self,para_path="/home/zijinping/Desktop/zijinping/resources/wy.para"):
        self.dict={}
        with open(para_path) as f:
            lines = f.readlines()
            for line in lines:
                line = line.rstrip()
                if len(line)==0 or line[0]=='#' or line[:4]=="gmt ": # ignore comment line and gmt set line
                    continue
                if line[:9]=="root_path":
                     self.dict["root_path"] = re.split("=",line.rstrip())[1]
                     continue
                content = re.split(" +",line.rstrip())[0]
                para,info = re.split("=",content)
                if len(re.split("\$",info))>1: # $ indicates citation of other parameters
                    for seg in re.split("\$",info)[1:]: # seperate each citation
                        sub = re.split("[/]",seg)[0]   # get the cited parameter name
                        info=info.replace("$"+sub,self.dict[sub])
                self.dict[para]=info
        f.close()
        #===============================================================================================
        # load in information
        # load_list=["city_loc",'city_label','ml_fault','zg_fault','Neo_fault','sta_loc','well']
        tmp_arr = []
        with open(self.dict['ml_fault'],'r') as f:
            for line in f:
                line = line.rstrip()
                _lon,_lat = re.split(" +",line)
                tmp_arr.append([float(_lon),float(_lat)])
        f.close()
        self.dict['ml_fault']=tmp_arr
        #------------------------------------------------------------------------------------------------
        tmp_dict={}
        count = 0
        with open(self.dict['zg_faults'],'r') as f:
            for line in f:
                line = line.rstrip()
                if line[0] == "#":
                    continue     # pass comment line
                elif line[0]==">":
                    count+=1
                    tmp_dict[count]=[]
                else:
                    _lon,_lat = re.split(" +",line)
                    tmp_dict[count].append([float(_lon),float(_lat)])
        f.close()
        self.dict['zg_faults'] = tmp_dict
        #---------------------------------------------------------
        tmp_dict={}
        count = 0
        with open(self.dict['Neo_faults'],'r') as f:
            for line in f:
                line = line.rstrip()
                if line[0] == "#":
                    continue     # pass comment line
                elif line[0]==">":
                    count+=1
                    tmp_dict[count]=[]
                else:
                    _lon,_lat = re.split(" +",line)
                    tmp_dict[count].append([float(_lon),float(_lat)])
        f.close()
        self.dict['Neo_faults'] = tmp_dict
        #---------------------------------------------------------
        tmp_arr = []
        with open(self.dict['city_locs'],'r') as f:
            for line in f:
                line = line.rstrip()
                _lon,_lat,_lvl,name = re.split(" +",line)[:4]
                tmp_arr.append([float(_lon),float(_lat),int(_lvl),name])
        f.close()
        self.dict['city_locs']=tmp_arr
        f.close()
        #----------------------------------------------------------------
        tmp_arr = []
        with open(self.dict['sta_locs'],'r') as f:
            for line in f:
                line = line.rstrip()
                net,sta,_lon,_lat,_ele,marker = re.split(" +",line)
                tmp_arr.append([float(_lon),float(_lat),float(_ele),net,sta,marker])
        f.close()
        self.dict["sta_locs"]=tmp_arr
        f.close()
        #---------------------------------------------------------------
        tmp_arr = []
        with open(self.dict['wells'],'r') as f:
            for line in f:
                line = line.rstrip()
                _lon,_lat,name,marker = re.split(" +",line)
                tmp_arr.append([float(_lon),float(_lat),name,marker])
        f.close()
        self.dict["wells"]=tmp_arr
        f.close()

def sta2inv(sta_file,out_file):
    net_stas = []
    cont = []
    with open(sta_file,'r') as f:
        for line in f:
            line = line.rstrip()
            net,sta,_lon,_lat,_ele,label = re.split(" +",line)
            net_sta = net+"_"+sta
            net_stas.append(net_sta)
            cont.append([net,sta,float(_lon),float(_lat),int(_ele),label])
    f.close()
    f_inv = open(out_file,'w')

    for tmp in cont:
        lat = tmp[3]
        lon = tmp[2]
        ele = tmp[4]
        net = tmp[0]
        sta = tmp[1]
        label = tmp[5]
        net_sta = net+sta
        lon_i = int(lon)
        lon_f = lon-lon_i
        lat_i = int(lat)
        lat_f = lat-lat_i
        f_inv.write(format(sta,"<6s")+format(net,"<4s")+"SHZ  "+format(lat_i,">2d")+" "+\
                    format(lat_f*60,">7.4f")+" "+format(lon_i,">3d")+" "+format(lon_f*60,">7.4f")+\
                    "E"+format(ele,">4d")+"\n")
    f_inv.close()


def sta2dd(sta_file,out_file):
    net_stas = []
    cont = []
    with open(sta_file,'r') as f:
        for line in f:
            line = line.rstrip()
            net,sta,_lon,_lat,_ele,label = re.split(" +",line)
            net_sta = net+"_"+sta
            net_stas.append(net_sta)
            cont.append([net,sta,float(_lon),float(_lat),int(_ele),label])
    f.close()
    f_dd = open(out_file,'w')

    for tmp in cont:
        lat = tmp[3]
        lon = tmp[2]
        ele = tmp[4]
        net = tmp[0]
        sta = tmp[1]
        label = tmp[5]
        net_sta = net+sta
        lon_i = int(lon)
        lon_f = lon-lon_i
        lat_i = int(lat)
        lat_f = lat-lat_i
        f_dd.write(format(net_sta,"<9s")+format(lat_i+lat_f,">9.6f")+format(lon_i+lon_f,">12.6f")+\
                   " "+format(ele,'>5d')+"\n")
    f_dd.close()
