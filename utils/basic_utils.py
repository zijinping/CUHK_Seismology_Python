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

def pha_subset(pha_file,loc_filter,obs_filter=8):
    """
    *.pha file is the input file for hypoDD ph2dt, this function subset the
    pha file by the boundary condition and the minimum observation condition.
    The output file is a file with ".st" suffix

    Parameters
    ----------
    pha_file: Str. The input file.
    loc_filter: array in format [lon_min, lon_max, lat_min, lat_max]
    obs_filter: The minimum observation
    """

    lon_min, lon_max, lat_min, lat_max = loc_filter
    out_file = pha_file+".st"
    f = open(out_file,"w")
    f.close()
    pha_content = []
    with open(pha_file,"r") as f:
        for line in f:
            pha_content.append(line.rstrip())
    f.close()
    i = 0
    j = 0
    record_list=[]
    for line in pha_content:
        if line[0]=="#":
            if i>0 and len(record_list) > (obs_filter+1):
                j=j+1
                with open(out_file,"a") as f:
                    for record in record_list:
                        f.write(record+"\n")
                f.close()
                record_list = []
                record_list.append(line)
            else:
                record_list = []
                record_list.append(line)
            i=i+1
            lat = float(re.split(" +",line)[7])
            lon = float(re.split(" +",line)[8])
            if lat>lat_min and lat<lat_max and lon>lon_min and lon<lon_max:
                region_pass = True
            else:
                region_pass = False
        else:
            if region_pass:
                record_list.append(line)
    if i>0 and len(record_list) > (obs_filter+1):
        j=j+1
        with open(out_file,"a") as f:
            for record in record_list:
                f.write(record+"\n")
        f.close()
    print("Event before filtering",i)
    print("Events qty after filtering",j)
     

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
        if t2 < starttime or t1 > endtime:
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
    sec = string[12:]
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

def load_sum(sum_file):
    """
    *.sum file is the catalog summary file after Hyperinverse.
    This function returns a dict:
        -key is event id
        -value is an array with below component:
            --Str format event time "yyyymmddhhmmss**", also the event folder.
            --event longitude
            --event latitude
            --event depth
            --event magnitude
    """
    sum_list = {}
    with open(sum_file,'r') as f:
        for line in f:
            eve_id=int(line[136:146])
            eve_folder = line[0:16]
            evla = int(line[16:18])+0.01*int(line[19:23])/60
            evlo = int(line[23:26])+0.01*int(line[27:32])/60
            evdp = int(line[31:36])*0.01
            e_mag = int(line[123:126])*0.01
            sum_list[eve_id] = [eve_folder,evlo,evla,evdp,e_mag]
    return sum_list

def load_sum_rev(sum_file):
    """
    *.sum file is the catalog summary file after Hyperinverse.
    This function returns a dict:
        -key is event time in "yyyymmddhhmmss**" format, same with event folder
        -value is an array with below component:
            --event id
            --event longitude
            --event latitude
            --event depth
            --event magnitude
    """

    sum_list = {}
    with open(sum_file,'r') as f:
        for line in f:
            eve_id=int(line[136:146])
            eve_folder = line[0:16]
            evla = int(line[16:18])+0.01*int(line[19:23])/60
            evlo = int(line[23:26])+0.01*int(line[27:32])/60
            evdp = int(line[31:36])*0.01
            e_mag = int(line[123:126])*0.01
            sum_list[eve_folder] = [eve_id,evlo,evla,evdp,e_mag]
    f.close()
    return sum_list

def load_hypoDD(reloc_file="hypoDD.reloc",shift_hour=0):
    """
    load results of hypoDD
    return eve_dict, df

    Parameters
    ----------
    If the time of results is not in UTC time zone, a time shift is needed.
    For example, Beijing time zone is 8 hours early than UTC time, 8 hours 
    should be deducted so as to be consistent with UTC time.
    """

    eve_dict={}
    columns = ["ID","LAT","LON","DEPTH","X","Y","Z","EX","EY","EZ",\
           "YR","MO","DY","HR","MI","SC","MAG",\
           "NCCP","NCCS","NCTP","NCTS","RCC","RCT","CID"]
    number = 0
    with open(reloc_file,"r") as f:
        for line in f:
            number = number+1
            if number%10==0:
                print("Current in process %d    "%number,end='\r')
            data=re.split(" +",line.rstrip())[1:]
            try:
                data_arr = np.vstack((data_arr,data))
            except:
                data_arr = np.array(data)
            eve_id = data[0]
            eve_lat = data[1]
            eve_lon = data[2]
            eve_dep = data[3]
            year = int(data[10])
            month = int(data[11])
            day = int(data[12])
            hour = int(data[13])
            minute = int(data[14])
            seconds = float(data[15])
            eve_time = UTCDateTime(year,month,day,hour,minute)+float(data[15])-shift_hour*60*60
            eve_time_str=str(eve_time)
            eve_mag =data[16]
            eve_dict[eve_time_str]=[float(eve_lon),float(eve_lat),float(eve_dep),float(eve_mag),int(eve_id)]
    f.close()
    df = pd.DataFrame(data=data_arr,columns=columns)
    return eve_dict,df

def hypoDD_mag_mapper(reloc_file,out_sum):
    """
    The output of hypoDD doesn't contain magnitude information.
    This function reads magnitude information from *.sum file, which is the
    output of hyperinverse and provide to hypoDD file.
    
    The results will cover the input reloc_fiie
    """

    #get the magnitude dictionary
    event_mag_list = {}
    with open(out_sum,"r") as f_obj:
        for line in f_obj:
            event_id = int(line[136:146])
            event_mag = int(line[122:126])*0.01
            event_mag_list[event_id]=event_mag
    f_obj.close()
    #add in the magnitude
    new_dd = []
    with open(reloc_file,"r") as f_obj:
        for line in f_obj:
            dd_event_id = int(line[0:9])
            dd_event_mag = event_mag_list[dd_event_id]
            new_line=line[:128]+format(dd_event_mag,'5.2f')+line[132:]
            new_dd.append(new_line)
    f_obj.close()
    with open(reloc_file,"w") as f_obj:
        for line in new_dd:
            f_obj.write(line)
    f_obj.close()

def hypoDD_ref_days(reloc_file,ref_time,shift_hours=0):
    """
    Add one column to the last of hypoDD files, calculate the length of time 
    between the referece time and the event time in days.
    The output is a file with the same title with reloc_file and add ".add" as
    suffix.

    Parameters
    ----------
     reloc_file: The hypoDD relocation file.
       ref_time: Reference time in UTCDateTime format
    shift_hours: Used when event time is not in UTC time zone
    """

    new_add=[]
    with open(reloc_file,"r") as f:
        for line in f:
            year = int(re.split(" +",line)[11])
            month = int(re.split(" +",line)[12])
            day = int(re.split(" +",line)[13])
            hour = int(re.split(" +",line)[14])
            minute = int(re.split(" +",line)[15])
            seconds = float(re.split(" +",line)[16])
            eve_time = UTCDateTime(year,month,day,hour,minute,0)+seconds - shift_hours*60*60
            days = (eve_time - ref_time)*1.0/(24*60*60)
            new_line=line[:-1]+" "+format(days,'4.2f')
            new_add.append(new_line)
    f.close()
    with open(reloc_file+".add","w") as f:
        for line in new_add:
            f.write(line+"\n")
    f.close()

