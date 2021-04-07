#/usr/bin/bash
###################################
# coding: utf-8
# Author: ZI,Jinping
# History:
#     2021-04-06 Initial coding
###################################

import os
import numpy as np
from math import radians,cos,acos,sin,asin,sqrt,ceil,pi,floor
import obspy
from obspy import Stream
import glob
import re
from obspy import UTCDateTime
import pandas as pd
from utils.basic_utils import str2time,time2str
import matplotlib.pyplot as plt

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
            eve_time = UTCDateTime(year,month,day,hour,minute)+seconds-shift_hour*60*60
            eve_time_str=time2str(eve_time)
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

def hypoDD_hist(dd_file="hypoDD.reloc",ref_time=UTCDateTime(2019,3,1,0,0,0)):
    """
    Plot hypoDD results in a histogram plot.
    
    Parameters:
    dd_file: Path of hypoDD file
    ref_time: Reference time for plot
    """
    eve_list = list(eve_dict)
    ref_list = []
    time_list = []
    number=0
    with open(dd_file,"r") as f:
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
            eve_time = UTCDateTime(year,month,day,hour,minute)+seconds
            time_list.append(eve_time)
            gap_time = eve_time - ref_time
            ref_list.append((eve_time-ref_time)/(60*60*24))
    min_day=floor(min(ref_list))
    max_day=ceil(max(ref_list))
    bins = np.linspace(min_day,max_day,max_day-min_day+1)
    fig1 = plt.figure(1,figsize=(8,4))
    ax1 = plt.subplot(1,1,1)
    ax1.hist(ref_list,bins)
    # The bottom x-axis is in days
    ax1.set_xlim([0,max_day])
    # The top x-axis marks year and month in YYYYMM
    tick_list_1 = [] # Store the position number
    tick_list_2 = [] # Store the tick text
    ref_year = ref_time.year
    ref_month = ref_time.month
    ref_day = ref_time.day
    if ref_day == 1:
        tick_list_1.append(0)
        tick_list_2.append(str(ref_year)+str(ref_month).zfill(2))
    status = True # Start to loop month by month
    loop_time = UTCDateTime(ref_year,ref_month,1) # Initiate loop time
    step = 32 #32 > 31. Make sure each step pass to next month
    while status==True:
        loop_time = loop_time + step*24*60*60
        tmp_year = loop_time.year
        tmp_month = loop_time.month
        loop_time = UTCDateTime(tmp_year,tmp_month,1)
        diff_days = (loop_time - ref_time)/(24*60*60)
        if diff_days > (max_day):
            status=False
        else:
            tick_list_1.append(diff_days)
            tick_list_2.append((str(tmp_month).zfill(2)))
    ax2 = ax1.twiny()
    ax2.set_xlim([0,max_day])
    ax2.plot(0,0,'k.')
    plt.xticks(tick_list_1,tick_list_2)
    ax1.set_xlabel("Time, days")
    ax1.set_ylabel("event quantity")
    ax2.set_xlabel("date")
    plt.show()
