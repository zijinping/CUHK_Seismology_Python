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

def str2time(string):
    """
    Convert from "yyyymmddhhmmss**" to UTCDateTime and return
    """
    year = string[0:4]
    mo = string[4:6]
    dy = string[6:8]
    hr = string[8:10]
    min = string[10:12]
    sec = string[12:18]
    time = UTCDateTime(year+"-"+mo+"-"+dy+"T"+hr+":"+min+":"+sec[0:2]+"."+sec[2:])
    return time

def time2str(time):
    """
    Convert from UTCDateTime to "yyyymmddhhmmss**" format and return

    Parameter
    ---------
    time: UTCDateTime format
    """

    year = time.year
    month = time.month
    day = time.day
    hour = time.hour
    minute = time.minute
    second = time.second
    m_second = time.microsecond
    string = str(year)+str(month).zfill(2)+str(day).zfill(2)+str(hour).zfill(2)+\
            str(minute).zfill(2)+str(second).zfill(2)+str(m_second).zfill(6)
    return string

